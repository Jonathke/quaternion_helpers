from helpers import *


if __name__ == "__main__":
    p = 10007
    B = QuaternionAlgebra(-1, -p)
    O0 = B.maximal_order()

    #Generate random orders, compute their successive minima
    for _ in range(10):
        I = heuristic_random_ideal(O0, next_prime(randint(p, p**2)))
        O = I.right_order()

        mins = successive_minima(O)
        print("\n~~~~~~~~~~~~~~~~")
        print(O)
        print(f"successive minima: {mins}")
        print(f"successive minima log_p: {[round(log(lam, p), 5) for lam in mins]}")
