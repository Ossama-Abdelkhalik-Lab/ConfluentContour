"""Multi-revolution tables: T_min(lam, m), v*(lam, m) and the branch roots
v(lam, sigma; m), sigma = branch * sqrt(1 - (T_min/t)^(1/3)) in [-1, 1]."""
import numpy as np
S2 = np.sqrt(2.0); PI = np.pi; TM = 0.995

def geom(lam):
    lam = np.asarray(lam, float)
    d = np.where(lam >= 0, 1.0, -1.0)
    sh = (1 - lam * lam) / (1 + lam * lam); cs = np.sqrt(np.maximum(0, 1 - sh * sh)); ss = sh
    a = d * S2 * cs; b = S2 * a
    y0 = np.where(d > 0, 2 * ss * ss / (1 + cs), 2 * (1 + cs)); yend = 2 * ss * ss / (1 + cs)
    return a, b, y0, yend

def T_m(v, m, a, b, y0, yend):
    sh, ch = np.sin(v / 2), np.cos(v / 2)
    y = np.where(b < 0, yend - 2 * b * ch * ch, y0 + 2 * b * sh * sh)
    sv, cv = np.sin(v), np.cos(v)
    return np.sqrt(y) * (y * (v + m * PI - sv * cv) / (S2 * sv ** 3) + a)

def dT(v, m, a, b, y0, yend, h=1e-6):
    return (T_m(v + h, m, a, b, y0, yend) - T_m(v - h, m, a, b, y0, yend)) / (2 * h)

def tmin(lam, m):
    a, b, y0, yend = geom(lam)
    lo = np.full(np.shape(lam), 1e-3); hi = np.full(np.shape(lam), PI - 1e-3)
    for _ in range(60):                       # bisection on dT/dv (unimodal T)
        mid = .5 * (lo + hi); f = dT(mid, m, a, b, y0, yend) < 0
        lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
    vs = .5 * (lo + hi)
    return T_m(vs, m, a, b, y0, yend), vs

def branch_root(lam, sig, m):
    a, b, y0, yend = geom(lam); tm, vs = tmin(lam, m)
    t = tm / (1 - sig * sig) ** 3
    lower = sig < 0
    lo = np.where(lower, 1e-12, vs); hi = np.where(lower, vs, PI - 1e-12)
    for _ in range(80):
        mid = .5 * (lo + hi); Tv = T_m(mid, m, a, b, y0, yend)
        right = np.where(lower, Tv > t, Tv < t)          # root lies right of mid
        lo = np.where(right, mid, lo); hi = np.where(right, hi, mid)
    return .5 * (lo + hi), vs

def lag4(x, x0, h, n):
    u = (x - x0) / h; i = np.clip(np.floor(u).astype(int) - 1, 0, n - 4); p = u - i
    return i, np.stack([-(p-1)*(p-2)*(p-3)/6, p*(p-2)*(p-3)/2, -p*(p-1)*(p-3)/2, p*(p-1)*(p-2)/6])

def interp2(Tab, n1, n2, x, y, x0, x1, y0, y1):
    i, wi = lag4(x, x0, (x1 - x0) / (n1 - 1), n1); j, wj = lag4(y, y0, (y1 - y0) / (n2 - 1), n2)
    out = 0
    for a in range(4):
        for b in range(4):
            out = out + wi[a] * wj[b] * Tab[i + a, j + b]
    return out

if __name__ == "__main__":
    rng = np.random.default_rng(3)
    N = 20000
    r2 = rng.uniform(0.1, 1, N); th = rng.uniform(1e-3, PI - 1e-3, N); d = rng.choice([-1, 1], N)
    c = np.sqrt((1 - r2) ** 2 + 4 * r2 * np.sin(th / 2) ** 2); sl = (1 + r2 + c) / 2
    lam = d * np.sqrt(1 - c / sl); tau = lam / (1 + np.sqrt(c / sl))
    sig = rng.uniform(-0.999, 0.999, N)
    for m in (1, 2, 5):
        vt, vs_t = branch_root(lam, sig, m); tm_t, _ = tmin(lam, m)
        dist = np.minimum(vt, PI - vt)
        for (nl, ns) in ((32, 32), (48, 48), (64, 64)):
            Tg = np.linspace(-TM, TM, nl); Lg = 2 * Tg / (1 + Tg * Tg); Sg = np.linspace(-1, 1, ns)
            LL, SS = np.meshgrid(Lg, np.clip(Sg, -1 + 1e-9, 1 - 1e-9), indexing='ij')
            Vtab, _ = branch_root(LL, SS, m)
            tmg, vsg = tmin(Lg, m)
            ve = interp2(Vtab, nl, ns, tau, sig, -TM, TM, -1, 1)
            i, w = lag4(tau, -TM, 2 * TM / (nl - 1), nl)
            tme = sum(w[k] * tmg[i + k] for k in range(4))
            rel = np.abs(ve - vt) / dist; rt = np.abs(tme / tm_t - 1)
            kb = (nl * ns + 2 * nl) * 8 / 1024
            print(f"m={m} grid {nl}x{ns} ({kb:5.1f} KB): v rel med {np.median(rel):.1e} p99 {np.quantile(rel,.99):.1e} "
                  f"frac<3e-5 {np.mean(rel<3e-5):.3f} | Tmin rel max {rt.max():.1e} p99 {np.quantile(rt,.99):.1e}")
