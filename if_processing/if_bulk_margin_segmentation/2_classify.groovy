// ============================================================================
// Script 2: tumour-cell classification  (transform -> Otsu | GMM)
// ----------------------------------------------------------------------------
// Splits cells into Tumour+ / Tumour- on one intensity measurement.
//
// Two knobs matter, in this order:
//   transform : reshapes the distribution. Fluorescence is right-skewed and ~log-
//               normal per population, so 'log' or 'arcsinh' separates the two humps
//               and makes the threshold well-defined. THIS is what changes the split.
//   method    : 'otsu' (nonparametric) or 'gmm' (fit gmmComponents Gaussians;
//               ONLY the highest-mean component is Tumour, cut where it overtakes
//               the 2nd-highest). GMM suits log/asinh-transformed IF data.
//
// NOTE: linear scaling (z-score / min-max) does NOT change the result on one image -
// Otsu and GMM both partition at a boundary that scales with the data, so an affine
// transform moves values and boundary together. Use z-scoring only to share ONE
// cutoff across a batch (set method='fixed', cut=<value in transformed space>).
// ============================================================================

// ---- EDIT: exact measurement name (check one cell in the Measurements table) ----
def measName  = 'ROI: 1.00 \u00B5m per pixel: tdtomato: Mean' //<---- check the name of the measurement to use here.
//def measName = 'Cell: tdtomato: Mean'

def transform = 'log'      // 'none' | 'log' | 'arcsinh'
def method    = 'gmm'      // 'otsu' | 'gmm' | 'fixed'
double cofactor = 1.0      // arcsinh cofactor ~ noise scale (ignored for log/none)
double fixedCut = Double.NaN // used only when method='fixed' (value in TRANSFORMED space)
// ---------------------------------------------------------------------------

def cells = getDetectionObjects().findAll { it.getMeasurementList().containsKey(measName) }
def raw = cells.collect { it.getMeasurementList().get(measName) }
                .findAll { it != null && Double.isFinite(it) } as double[]
if (raw.length == 0) { println "No values for '${measName}' - check the name."; return }

// --- transform --------------------------------------------------------------
// log: shift so the smallest value maps to ~0, then log1p (handles background-
// subtracted zeros/negatives gracefully). arcsinh: asinh(v/cofactor), also fine
// with negatives. The SAME transform is applied to each cell below when classifying.
double gmin = raw.min()
def tf = { double v ->
    if (transform == 'log')     return Math.log1p(Math.max(0.0, v - gmin))
    if (transform == 'arcsinh') return Math.log(v/cofactor + Math.sqrt((v/cofactor)*(v/cofactor) + 1.0))
    return v
}
double[] values = raw.collect { tf(it) } as double[]

// --- Otsu (on transformed values) ------------------------------------------
def otsu = { double[] v ->
    double mn = v.min(); double mx = v.max()
    int nbins = 256
    double bw = (mx - mn) / nbins
    if (bw <= 0) return mn
    int[] hist = new int[nbins]
    v.each { int b = Math.min(nbins-1, (int)((it - mn)/bw)); hist[b]++ }
    double total = v.length, sum = 0
    for (int i=0;i<nbins;i++) sum += (mn + (i+0.5)*bw) * hist[i]
    double sumB=0, wB=0, best=0, thr=mn
    for (int i=0;i<nbins;i++){
        wB += hist[i]; if (wB==0) continue
        double wF = total - wB; if (wF==0) break
        double mB = sumB/wB
        sumB += (mn + (i+0.5)*bw) * hist[i]
        double mF = (sum - sumB)/wF
        double between = wB*wF*(mB-mF)*(mB-mF)
        if (between > best){ best = between; thr = mn + (i+1)*bw }
    }
    return thr
}

// --- K-component 1-D Gaussian mixture, fit by EM ----------------------------
// Fits gmmComponents Gaussians. ONLY the highest-mean component is Tumour, so the
// cutoff is the boundary where the top component overtakes the second-highest
// (i.e. the x above which the top component is the most likely). With K=3 this
// separates e.g. background / stroma / tumour and calls only the brightest tumour.
int gmmComponents = 4   // number of Gaussians to fit (>=2). Top one = Tumour.
def gauss = { double x, double m, double s -> Math.exp(-((x-m)*(x-m))/(2*s)) / Math.sqrt(2*Math.PI*s) }
def gmm = { double[] v ->
    def sorted = (v as List).sort()
    int n = sorted.size()
    int K = Math.max(2, gmmComponents)
    if (n < K) return sorted[(int)(n/2)]
    // init: means at K evenly spaced quantiles, equal weights, shared spread
    double[] m = new double[K], s = new double[K], w = new double[K]
    double gmean = sorted.sum()/n
    double gvar  = 0; sorted.each { gvar += (it-gmean)*(it-gmean) }; gvar = Math.max(gvar/n, 1e-9)
    for (int k=0;k<K;k++) {
        m[k] = sorted[(int)Math.min(n-1, (long)((k+0.5)/K * n))]
        s[k] = gvar; w[k] = 1.0/K
    }
    double[][] r = new double[K][n]
    for (int iter=0; iter<200; iter++) {
        // E-step: responsibilities
        for (int i=0;i<n;i++) {
            double x = sorted[i]; double den = 0
            for (int k=0;k<K;k++){ r[k][i] = w[k]*gauss(x,m[k],s[k]); den += r[k][i] }
            if (den<=0) den = 1e-300
            for (int k=0;k<K;k++) r[k][i] /= den
        }
        // M-step
        for (int k=0;k<K;k++) {
            double sr=0, sm=0
            for (int i=0;i<n;i++){ sr += r[k][i]; sm += r[k][i]*sorted[i] }
            double nm = sm/Math.max(sr,1e-9), ss=0
            for (int i=0;i<n;i++){ double d=sorted[i]-nm; ss += r[k][i]*d*d }
            m[k]=nm; s[k]=Math.max(ss/Math.max(sr,1e-9),1e-9); w[k]=sr/n
        }
    }
    // sort components by mean (ascending) so index K-1 is the Tumour component
    def order = (0..<K).toList().sort(false) { m[it] }   // toList(): a range is immutable, in-place sort throws
    double[] mo=new double[K], so=new double[K], wo=new double[K]
    for (int k=0;k<K;k++){ mo[k]=m[order[k]]; so[k]=s[order[k]]; wo[k]=w[order[k]] }
    def sb = new StringBuilder("  GMM (" + K + " comp):")
    for (int k=0;k<K;k++) sb << String.format(" [%d] N(%.3f, %.3f) w=%.2f", k, mo[k], Math.sqrt(so[k]), wo[k])
    println sb.toString()
    // cut = boundary where the top component (K-1) overtakes the 2nd-highest (K-2),
    // found by bisection between their means. Above this x the top comp is argmax.
    double loB=mo[K-2], hiB=mo[K-1]
    for (int it2=0; it2<100; it2++) {
        double x=(loB+hiB)/2
        double pTop = wo[K-1]*gauss(x,mo[K-1],so[K-1])
        double pSec = wo[K-2]*gauss(x,mo[K-2],so[K-2])
        if (pTop >= pSec) hiB = x else loB = x
    }
    return (loB+hiB)/2
}

// --- pick the cutoff --------------------------------------------------------
double cut
if      (method == 'fixed') { cut = fixedCut }
else if (method == 'gmm')   { cut = gmm(values) }
else                        { cut = otsu(values) }
println String.format("Method=%s transform=%s | cutoff (transformed) = %.4f  [raw min=%.3f]",
                      method, transform, cut, gmin)

// --- classify ---------------------------------------------------------------
def cTum = getPathClass('Tumour')
def cNeg = getPathClass('Stroma')
cells.each {
    def v = it.getMeasurementList().get(measName)
    boolean pos = (v != null && Double.isFinite(v) && tf(v) >= cut)
    it.setPathClass(pos ? cTum : cNeg)
}
fireHierarchyUpdate()
println "Tumour+ cells: " + cells.count { it.getPathClass() == cTum }
println "Tumour- cells: " + cells.count { it.getPathClass() == cNeg }