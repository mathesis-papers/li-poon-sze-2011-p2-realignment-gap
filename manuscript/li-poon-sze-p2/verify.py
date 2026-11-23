import json
from itertools import combinations
from math import comb, isqrt
from fractions import Fraction

import sympy as symbolic


def equal(left, right):
    assert symbolic.simplify(left - right) == 0, (left, right)


def square_root_interval(value, digits=70):
    scale = 10 ** digits
    value = Fraction(value)
    floor = isqrt(value.numerator * scale * scale // value.denominator)
    return Fraction(floor, scale), Fraction(floor + 1, scale)


def leading_interval(m, n, family):
    if family == 'main':
        lower, upper = square_root_interval(m * n)
        return 1 / upper, 1 / lower
    if family == 'supported':
        lower, upper = square_root_interval(m * m - 1)
        return 1 / (m * upper), 1 / (m * lower)
    return Fraction(1, m * m), Fraction(1, m * m)


def profile_interval(m, degree, interval, tail_interval=None):
    lower, upper = interval
    if tail_interval is None:
        tail_lower = (1 - upper) / (m * m - 1)
        tail_upper = (1 - lower) / (m * m - 1)
    else:
        tail_lower, tail_upper = tail_interval
    coefficient = comb(m * m - 1, degree) if degree < m * m else 0
    other = comb(m * m - 1, degree - 1)
    return (coefficient * tail_lower ** degree + lower * other * tail_lower ** (degree - 1),
            coefficient * tail_upper ** degree + upper * other * tail_upper ** (degree - 1))


def witness(m, n, family, noise=Fraction(0)):
    if family == 'main':
        sectors = n // m
        leading = 1 / symbolic.sqrt(m * n)
    elif family == 'supported':
        sectors = m * m - 1
        leading = 1 / (m * symbolic.sqrt(sectors))
    else:
        sectors = m * m
        leading = symbolic.Rational(1, m * m)
    remainder = n - m * sectors
    repeated = (1 - leading) / (m * m - 1)
    entangled = symbolic.simplify(repeated / symbolic.sqrt(sectors))
    occupied = symbolic.simplify(leading ** 2 - repeated / (m * symbolic.sqrt(sectors)))
    unused = leading ** 2 if family == 'main' else symbolic.S.Zero
    assert min(remainder, sectors - 1) >= 0
    for coefficient in [entangled, occupied, unused]:
        assert coefficient.is_nonnegative is True
    state = symbolic.MutableSparseMatrix(m * n, m * n, {})
    for sector in range(sectors):
        support = [alice * n + alice * sectors + sector for alice in range(m)]
        for row in support:
            for column in support:
                state[row, column] += entangled
    for alice in range(m):
        for bob in range(n):
            index = alice * n + bob
            state[index, index] += occupied if bob < m * sectors else unused
    leading_bounds = leading_interval(m, n, family)
    repeated_bounds = ((1 - leading_bounds[1]) / (m * m - 1),
                       (1 - leading_bounds[0]) / (m * m - 1))
    if noise:
        weight = symbolic.Rational(noise.numerator, noise.denominator)
        state = (1 - weight) * state + weight / (m * n) * symbolic.eye(m * n)
        entangled = (1 - weight) * entangled
        occupied = (1 - weight) * occupied + weight / (m * n)
        unused = (1 - weight) * unused + weight / (m * n)
        leading = symbolic.sqrt((1 - weight) ** 2 * leading ** 2 +
                                (2 * weight - weight ** 2) / (m * n))
        repeated = (1 - weight) * repeated
        lower_square = (1 - noise) ** 2 * leading_bounds[0] ** 2 + (2 * noise - noise ** 2) / (m * n)
        upper_square = (1 - noise) ** 2 * leading_bounds[1] ** 2 + (2 * noise - noise ** 2) / (m * n)
        leading_bounds = (square_root_interval(lower_square)[0],
                          square_root_interval(upper_square)[1])
        repeated_bounds = tuple((1 - noise) * bound for bound in repeated_bounds)
    assert state == state.H
    equal(symbolic.trace(state), 1)
    for alice in range(m):
        for other_alice in range(m):
            marginal_entry = sum(state[alice * n + bob, other_alice * n + bob]
                                 for bob in range(n))
            equal(marginal_entry, symbolic.Rational(int(alice == other_alice), m))
    realigned = symbolic.MutableSparseMatrix(m * m, n * n, {})
    for (row, column), entry in state.todok().items():
        alice, bob = divmod(row, n)
        other_alice, other_bob = divmod(column, n)
        realigned[alice * m + other_alice, bob * n + other_bob] = entry
    gram = realigned * realigned.H
    direction = symbolic.Matrix([int(row == column) for row in range(m) for column in range(m)])
    mu = sectors * entangled ** 2
    nu = sectors * (2 * entangled * occupied + m * occupied ** 2) + remainder * unused ** 2
    target = mu * symbolic.eye(m * m) + nu * direction * direction.T
    for row in range(m * m):
        for column in range(m * m):
            equal(gram[row, column], target[row, column])
    equal(mu, repeated ** 2)
    equal(mu + m * nu, leading ** 2)
    if noise:
        assert leading_bounds[1] + (m * m - 1) * repeated_bounds[1] < 1
        assert noise / (m * n) > 0
    else:
        equal(leading + (m * m - 1) * repeated, 1)
    tail_lower, tail_upper = square_root_interval(Fraction((m - 1) * (n - 1), m * n))
    threshold_bounds = 1 - tail_upper, 1 - tail_lower
    assert leading_bounds[1] < threshold_bounds[0]
    for degree in range(2, m * m + 1):
        assert profile_interval(m, degree, leading_bounds, repeated_bounds)[0] > profile_interval(m, degree, threshold_bounds)[1]
    return {'m': m, 'n': n, 'family': family, 'q': sectors, 'r': remainder,
            'noise_weight': str(noise),
            'state_eigenvalue_lower_bound': str(noise / (m * n)),
            'trace': 1, 'row_gram': 'exact equality',
            'alice_marginal': 'identity divided by m',
            'trace_norm': 'strictly below one' if noise else 1,
            'strict_degrees': list(range(2, m * m + 1))}


def symbolic_identities():
    m, q, r, a, b = symbolic.symbols('m q r a b', positive=True)
    t = b / symbolic.sqrt(q)
    u = a ** 2 - b / (m * symbolic.sqrt(q))
    mu = q * t ** 2
    nu = q * (2 * t * u + m * u ** 2) + r * a ** 4
    equal(mu, b ** 2)
    equal(mu + m * nu, a ** 4 * m * (m * q + r))
    equal(m * q * (t + m * u) + m * r * a ** 2, a ** 2 * m * (m * q + r))
    equal(m ** 4 * (m ** 2 - 2) - (m ** 3 - m - 1) ** 2,
          2 * m ** 3 - m ** 2 - 2 * m - 1)
    noise = symbolic.symbols('noise', nonnegative=True)
    noisy_entangled = (1 - noise) * t
    noisy_occupied = (1 - noise) * u + noise / (m * (m * q + r))
    noisy_unused = (1 - noise) * a ** 2 + noise / (m * (m * q + r))
    equal(q * (noisy_entangled + m * noisy_occupied) ** 2 + m * r * noisy_unused ** 2,
          (1 - noise) ** 2 * a ** 4 * m * (m * q + r) +
          2 * noise * (1 - noise) * a ** 2 + noise ** 2 / (m * (m * q + r)))
    x = symbolic.symbols('x')
    for size in range(4, 18):
        tail = (1 - x) / (size - 1)
        for degree in range(2, size + 1):
            first = comb(size - 1, degree) if degree < size else 0
            value = first * tail ** degree + x * comb(size - 1, degree - 1) * tail ** (degree - 1)
            derivative = symbolic.Rational(degree - 1, size - 1) * comb(size - 1, degree - 1) * tail ** (degree - 2) * (tail - x)
            equal(symbolic.diff(value, x), derivative)
    count = 0
    for dimension in range(2, 51):
        for sectors in range(1, dimension * dimension - 1):
            assert (dimension * (sectors + 1) - 1) ** 2 <= dimension ** 4 * sectors
            count += 1
    return {'gram_and_trace_identities': 'exact symbolic identities',
            'profile_derivative': {'sizes': [4, 17], 'all_degrees': True},
            'dimension_inequality': {'m': [2, 50], 'pairs': count}}


def main():
    identities = symbolic_identities()
    cases = [(2, 3, 'main'), (2, 5, 'main'), (2, 6, 'supported'),
             (3, 26, 'supported'), (2, 8, 'uniform'), (3, 27, 'uniform')]
    examples = [witness(*case) for case in cases]
    noisy_examples = [witness(*case, noise=Fraction(1, 1000)) for case in cases]
    equal_profile = [symbolic.Rational(1, 2)] + [symbolic.Rational(1, 6)] * 3
    equal_values = [sum(symbolic.prod(equal_profile[index] for index in subset)
                    for subset in combinations(range(4), degree))
                    for degree in range(1, 5)]
    assert equal_values == [1, symbolic.Rational(1, 3), symbolic.Rational(5, 108), symbolic.Rational(1, 432)]
    print(json.dumps({'status': 'pass', 'arithmetic': 'exact symbolic and rational interval',
                      'interval_decimal_digits': 70,
                      'identities': identities, 'examples': examples,
                      'full_rank_ccnr_interior_examples': noisy_examples,
                      'two_qubit_equal_maxima': list(map(str, equal_values))}, indent=2))


if __name__ == '__main__':
    main()
