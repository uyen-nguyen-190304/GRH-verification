"""
discriminants.py

Functions to check if a discriminant is fundamental

Functions:
    - is_square_free(n): Check if an integer is free of any squared prime divisors
    - is_fundamental_discriminant(d): Check if a discriminant is fundamental of a real quadratic field
"""

def is_square_free(n: int) -> bool:
    """
    Purpose:
        Determine whether the absolute value of n is square-free
        A nonzero integer n is square-free if no prime p satisfies p^2 | n
    Input:
      n (int): The integer to test (can be negative)
    Return:
        bool: True if |n| is square-free (including n = 1), False otherwise
    """
    # Reduce n to its absolute value |n|
    n = abs(n)

    # 0 is not square-free; 1 is the boundary case, considered square-free
    if n in (0, 1):
        return n == 1   
    
    # If divisible by 4, it has factor 2^2 -> not square-free
    if n % 4 == 0:
        return False    
    
    # Check odd primes and p=2
    p = 2
    while p * p <= n:
        # If p^2 divides n, then n is not square-free
        if n % (p * p) == 0:
            return False   
        # Advance to next potential prime factor: 2 -> 3, then skip evens
        p +=1 if p == 2 else 2  

    # No square divisors found
    return True


def is_fundamental_discriminant(d: int) -> bool:
    """
    Purpose:
        Check if d is a fundamental discriminant
        A discriminant d of a quadratic number field must satisfy either of:
            (1) d ≡ 1 mod 4 and |d| is square-free, or
            (2) d ≡ 0 mod 4, let q = d/4; q ≡ 2 or 3 mod 4, and |q| is square-free
    Input:
        d (int): A given discriminant
    Return:
        bool: True if d is a fundamental discriminant, False otherwise
    """
    # Early check: 0 is not a fundamental discriminant
    if d == 0:
        return False
    
    # Case 1: d ≡ 1 mod 4
    if d % 4 == 1:
        # Must be square-free
        return is_square_free(d)

    # Case 2: d ≡ 0 mod 4
    if d % 4 == 0:
        q = d // 4
        # q must be square-free and q ≡ 2 or 3 (mod 4)
        return is_square_free(q) and q % 4 in (2, 3)
    
    # Otherwise, not a fundamental discriminant
    return False
