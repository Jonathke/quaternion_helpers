from sage.all import *
proof.all(False)

from qmat_iso_q import QMatModq

def mat_crt(res, mod):
    """
    Input:
    - res: a list of 2x2 matrices M_i defined over Zn_i
    - mod: the list of the n_i
    Output:
    - a matrix M such that M mod n_i == M_i
    """
    a_i = [ZZ(Mi[0, 0]) for Mi in res]
    b_i = [ZZ(Mi[0, 1]) for Mi in res]
    c_i = [ZZ(Mi[1, 0]) for Mi in res]
    d_i = [ZZ(Mi[1, 1]) for Mi in res]
    a, b, c, d = [crt(x, mod) for x in [a_i, b_i, c_i, d_i]]
    N = product(mod)
    return Matrix(Integers(N), 2, [a, b, c, d])

def quat_crt(res, mod):
    """
    Input:
    - res: a list of 4-tuples, each representing the coordinates of a
      quaternion theta_i in the fixed basis mod n_i
    - mod: the list of the n_i
    Output:
    - the 4-tuple of coordinates of the quaternion element theta such that
      theta mod n_i = theta_i
      """
    a_i = [ZZ(ri[0]) for ri in res]
    b_i = [ZZ(ri[1]) for ri in res]
    c_i = [ZZ(ri[2]) for ri in res]
    d_i = [ZZ(ri[3]) for ri in res]
    a, b, c, d = [crt(x, mod) for x in [a_i, b_i, c_i, d_i]]
    return (a, b, c, d)

class QMatIso:
    """
    Given a prime p = 3 mod 4 and a number n compute the isomorphism between
    B/nB and M_2(n) explicitly.
    n can be either a prime or provided as a list of its factors in sage format
    [(p_i, e_i)].
    """
    def __init__(self, p, n, B=None, O0=None):
        """
        Input:
        - p: base prime
        - n: modulo (either a prime or a list of factors and exponents)
        - B: optional, quaternion algebra ramified at -1 and -p
        - O0: optional, the order (currently only standard O0 supported)
        """
        self.p = p

        if not B:
            B = QuaternionAlgebra(-1, -p)
        self.B = B
        i, j, k = B.gens()
        self.B_basis = [1, i, (i+j)/2, (1+k)/2]

        if not O0:
            O0 = B.maximal_order(order_basis=(B(1), i, (i+j)/2, (1-k)/2))
        else:
            raise NotImplementedError("only O0 supported for now")

        self.O0 = O0

        if type(n) in [int, Integer]:
            assert is_prime(n)
            self.n_fact = [(n, 1)]
        else:
            self.n_fact = n

        self.iso_mod_q = [QMatModq(p, q_i, e_i) for q_i, e_i in self.n_fact]
        self.z_mod_q = [QMIq.Zn for QMIq in self.iso_mod_q]
        self.mods = [ZZ(q_i**e_i) for q_i, e_i in self.n_fact]

    def quat_to_mat(self, alpha):
        """
        Given a quaternion alpha return the image of alpha in M_2(Zn) under the
        isomorphism. We are working with basis (1, i, (i+j)/2, (1+k)/2)
        """
        if not alpha in self.O0:
            raise ValueError(f"{alpha} not in the order")

        # Recover the matrix mod q_i
        Mi = [QMIq.quat_to_mat(alpha) for QMIq in self.iso_mod_q]
        M = mat_crt(Mi, self.mods)

        return M

    def mat_to_quat(self, M):
        """
        Given a matrix M in M_2(Zn) return a lift of M in B.
        """
        # Recover the coordinates w.r.t. the basis
        elt_i = [QMIq.mat_to_quat(M) for QMIq in self.iso_mod_q]
        elt = quat_crt(elt_i, self.mods)
        alpha = sum(self.B_basis[i] * ZZ(elt[i]) for i in range(4))

        # alpha = self.B(quat_crt(alpha_i, self.mods))
        assert alpha in self.O0, "output not in the order??"
        return alpha

    def project(self, alpha):
        return self.quat_to_mat(alpha)

    def lift(self, M):
        return self.mat_to_quat(M)


if __name__ == "__main__":
    p = 2**10*3**10*5**2*7*11 - 1

    B = QuaternionAlgebra(-1, -p)
    i,j,k = B.gens()

    q = 2
    e = 10
    N = q**e*3**10
    ZN = Integers(N)
    O0 = B.maximal_order(order_basis=(B(1), i, (i+j)/2, (1-k)/2))

    QMI = QMatIso(p, [(q, e), (3, 10)], B)

    # Quaternion to matrices (mod n)
    for _ in range(1):
        while True:
            alpha = B.random_element()
            if gcd(N, alpha.denominator()) != 1:
                continue
            alpha = B([ZN(i) for i in alpha])
            break
        M = QMI.quat_to_mat(alpha)
        alpha1 = B(QMI.mat_to_quat(M))
        print(f'{(alpha - alpha1)/N = }')
        assert (alpha - alpha1)/N in O0

    # Matrices to quaternions

    # Linearity for matrices

    # Linearity for quaternions



