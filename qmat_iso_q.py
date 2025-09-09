from sage.all import *
proof.all(False)

def hensel_solve(b, q, e):
    """
    Use Hensel lifting to solve the equation x^2 + y^2 + y = b modulo q^e.
    """
    # TODO: actually use hensel, also maybe check derivatives
    Zn = Integers(q**e)
    c1 = None
    for _ in range(10000):
        c2 = Zn.random_element()
        y = -b - c2**2 + c2
        if y.is_square() and (c2 != 0 or y != 0):
            c1 = sqrt(y)
            break
    assert c1, "c1 not found - stop being lazy and do that properly"
    assert c1**2 + c2**2 - c2 == -b
    return c1, c2

class QMatModq:
    """
    Contstruct the isomorphism for prime powers
    """
    def __init__(self, p, q, e, B=None, O0=None):
        """
        Input:
        - p: base prime
        - q: base prime for the isomorphism (different from p)
        - e: exponent
        - B: optional, quaternion algebra ramified at -1 and -p
        - O0: optional, the maximal order; for now only O0 is supported
        """
        self.p = p
        self.q = q
        self.e = e

        print(f'Starting for {q}^{e}')

        if not B:
            B = QuaternionAlgebra(-1, -p)
        self.B = B

        i, j, k = B.gens()
        self.B_basis = [1, i, (i+j)/2, (1+k)/2]

        # TODO: support for arbitrary order
        if not O0:
            O0 = B.maximal_order(order_basis=(B(1), i, (i+j)/2, (1-k)/2))
        else:
            raise NotImplementedError("arbitraty order not supported")

        self.O0 = O0

        n = q**e
        Zn = Integers(n)

        self.n = n
        self.Zn = Zn

        # Find 4 matrices representing 1, i, (i+j)/2, (1+k)/2
        One = Matrix(Zn, 2, [1, 0, 0, 1])
        I = Matrix(Zn, 2, [0, 1, -1, 0])
        assert I**2 == -One

        # (i+j)/2 has norm (p+1)/4 and trace zero
        b = Zn((p+1)/4)

        # Find c1, c2 s.t. c1^2 + c2^2 - c2 = b mod n
        c1, c2 = hensel_solve(b, q, e)

        J2 = Matrix(Zn, 2, [c1, c2, c2-1, -c1])

        assert J2**2 == -((p+1)//4) * One
        assert J2*I + I*J2 + One == 0

        # Now (1+k)/2 = 1 + i*(i+j)/2
        K2 = One + I * J2

        self.M_basis = [One, I, J2, K2]

    def quat_to_mat(self, alpha):
        """
        Given a quaternion alpha return the image of alpha in M_2(Zn) under the
        isomorphism. We are working with basis (1, i, (i+j)/2, (1+k)/2)
        """
        if not alpha in self.O0:
            raise ValueError(f"{alpha} not in the order")

        # Inverse of the basis matrix
        M = Matrix(QQ, 4, [
            [1, 0, 0, -1],
            [0, 1, -1, 0],
            [0, 0, 2, 0],
            [0, 0, 0, 2]])
        beta = M * vector(QQ, alpha)
        mat = sum(self.Zn(beta[i]) * self.M_basis[i] for i in range(4))
        return mat

    def mat_to_quat(self, M):
        """
        Given a matrix M in M_2(Zn) return a lift of M in B expressed in
        coordinates of the basis.
        """
        One, I, J2, K2 = self.M_basis
        c1, c2, c3, c4 = J2.list()
        k1, k2, k3, k4 = K2.list()
        A = Matrix(self.Zn, 4, [
            1, 0, c1, k1,
            0, 1, c2, k2,
            0, -1, c3, k3,
            1, 0, c4, k4
        ])
        v = Matrix(self.Zn, 4, 1, M.list())
        elt = A.solve_right(v).list()
        return elt
        # alpha = sum(self.B_basis[i] * ZZ(elt[i]) for i in range(4))
        # assert alpha in self.O0, "output not in the order??"
        # return alpha

    def project(self, alpha):
        return self.quat_to_mat(alpha)

    def lift(self, M):
        return self.mat_to_quat(M)

