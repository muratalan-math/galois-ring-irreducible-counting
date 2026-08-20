#include <bits/stdc++.h>
using namespace std;

// Directly enumerate all monic local polynomials over Z/8Z reducing to x^e,
// for 1 <= e <= 6. A polynomial is marked reducible iff it is a product of
// two positive-degree monic local factors. No manuscript recurrence or closed
// counting formula is used.

static vector<int> factor_from_code(uint32_t code, int deg) {
    vector<int> a(deg + 1, 0);
    for (int i = 0; i < deg; ++i) {
        a[i] = 2 * (code & 3u); // coefficients in 2(Z/8Z)={0,2,4,6}
        code >>= 2;
    }
    a[deg] = 1;
    return a;
}

static vector<int> multiply_mod8(const vector<int>& a, const vector<int>& b) {
    vector<int> c(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j)
            c[i + j] = (c[i + j] + a[i] * b[j]) & 7;
    return c;
}

static uint32_t encode_local(const vector<int>& c, int deg) {
    uint32_t code = 0, place = 1;
    for (int i = 0; i < deg; ++i) {
        int v = c[i] & 7;
        if (v & 1) {
            cerr << "internal error: coefficient not in 2(Z/8Z)\n";
            exit(3);
        }
        code += place * static_cast<uint32_t>((v / 2) & 3);
        place *= 4;
    }
    return code;
}

int main() {
    cout << "Direct local irreducibility over Z/8Z (residual x^e)\n";
    cout << "e  total  reducible  J_{2,3}(e)\n";
    for (int e = 1; e <= 6; ++e) {
        const uint32_t total = 1u << (2 * e); // 4^e local coefficient vectors
        vector<uint8_t> reducible(total, 0);

        for (int r = 1; r < e; ++r) {
            const uint32_t na = 1u << (2 * r);
            const uint32_t nb = 1u << (2 * (e - r));
            for (uint32_t ia = 0; ia < na; ++ia) {
                const auto a = factor_from_code(ia, r);
                for (uint32_t ib = 0; ib < nb; ++ib) {
                    const auto b = factor_from_code(ib, e - r);
                    reducible[encode_local(multiply_mod8(a, b), e)] = 1;
                }
            }
        }

        const uint32_t nred = accumulate(reducible.begin(), reducible.end(), 0u);
        cout << e << "  " << total << "  " << nred << "  " << (total - nred) << "\n";
    }
    return 0;
}
