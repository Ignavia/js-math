/*
    <EquationSolver: This program solves different kinds of equations.>
    Copyright (C) <2013>  <Lars Reimann>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

// =============================================================================
// Solvers
// =============================================================================

/**
 * Computes the set of all complex x that satisfy a = 0 and returns the result
 * as a string.
 *
 * @param  {Number} a
 * Corresponds to the a in a = 0.
 *
 * @return {String}
 * the set of all complex x that satisfy a = 0.
 */
function solveConstant(a) {
    return a === 0 ? "L = ℂ" : "L = ∅";
}

/**
 * Computes the set of all complex x that satisfy ax + b = 0 and returns the
 * result as a string.
 * 
 * @param  {Number} a 
 * Corresponds to the a in ax + b = 0.
 *
 * @param  {Number} b
 * Corresponds to the b in ax + b = 0.
 *
 * @param  {Boolean} numeric
 * If this parameter is true an approximate result is returned. Otherwise the
 * method tries to compute an exact result. This parameter is optional, the
 * default value is true. Note that computing an exact result can take
 * significantly more time.
 *
 * @return {String}
 * the set of all complex x that satisfy ax + b = 0.
 */
function solveLinear(a, b, numeric) {
    if (typeof numeric === "undefined") {
        numeric = true;
    }

    if (a === 0) {
        return solveConstant(b);
    } else if (numeric) {
        return "L = {" + (-b/a) + "}";
    } else {
        return "L = {" + new Fraction(-b, a) + "}";
    }
}

/**
 * Computes the set of all complex x that satisfy ax^2 + bx + c = 0 and prints
 * the result.
 *
 * @param  {Number} a
 * Corresponds to the a in ax^2 + bx + c = 0.
 * 
 * @param  {Number} b
 * Corresponds to the b in ax^2 + bx + c = 0.
 * 
 * @param  {Number} c
 * Corresponds to the c in ax^2 + bx + c = 0.
 * 
 * @param  {Boolean} numeric
 * If this parameter is true an approximate result is returned. Otherwise the
 * method tries to compute an exact result. This parameter is optional, the
 * default value is true. Note that computing an exact result can take
 * significantly more time.
 *
 * @return {String}
 * the set of all complex x that satisfy ax^2 + bx + c = 0.
 */
function solveQuadratic(a, b, c, numeric) {
    if (typeof numeric === "undefined") {
        numeric = true;
    }

    if (a === 0) {
        return solveLinear(b, c);
    } else if (numeric) {
        var disc = b*b - 4*a*c;
        var rat = -b/2*a;
        if (disc === 0) {
            return "L = {" + rat + "}";
        } else if (disc < 0) {
            return "L = {" + rat + " +/- " + abs(Math.sqrt(-disc)/(2*a)) + "*i}";
        } else {
            return "L = {" + rat + " +/- " + abs(Math.sqrt(disc)/(2*a)) + "}";
        }
    } else {
        var disc = b*b - 4*a*c;
        var rat = new Fraction(-b, 2*a);
        if (disc === 0) {
            return "L = {" + rat + "}";
        } else {
            return "L = {" + rat + " +/- " +
                   fracTimesRoot(new Fraction(1, 2*abs(a)), new Root(1, disc)) +
                   "}";
        }
    }
}

/**
 * Solves the system of linear equations. equation has to be a two dimensional
 * array of numbers. Each element of the array represents one equation that has
 * to be satisfied. Each equation is represented by the coefficients of the
 * unknown values. The leftmost value in each array is the coefficient of x0,
 * then comes the coefficient of x1 and so forth. The rightmost number is the
 * number in the righthandside of the equation. Example:
 * 
 * 3*x0 + 4*x1 - 1*x2 = 3 is represented as [3, 4, -1, 3].
 *
 * In order for this function to work, the system of equations must be fully
 * specified. Refer to the documentation of {@code isFullySpecified} for more
 * information. The result is returned as a string.
 *
 * @param  {Number[][]} equations
 * The system of equations. It has to be represented using the form specified
 * above.
 *
 * @param  {boolean} numeric
 * If this parameter is true an approximate result is returned. Otherwise the
 * method tries to compute an exact result. This parameter is optional, the
 * default value is true. Note that computing an exact result can take
 * significantly more time.
 *
 * @return {String}
 * A string representing the solution to this system of equations.
 */
function solveLinearSystem(equations, numeric) {
    if (typeof numeric === "undefined") {
        numeric = true;
    }

    if (isFullySpecified(equations)) {
        gaussianElimination(equations, numeric);
        var solution = solveLinearSystemHelper(equations, numeric);

        if (typeof solution === 'string') {
            return solution;
        }
        var s = "";
        for (var i = 0; i < solution.length; i++) {
            if (numeric) {
                s += "x" + i + " = " + linSysToStringHelperNumeric(solution[i]);
            } else {
                s += "x" + i + " = " + linSysToStringHelperExact(solution[i]);
            }
        }
        return s;
    } else {
        return "You did not specify all elements.";
    }
}

// =============================================================================
// Helper functions
// =============================================================================

/**
 * Returns a new null vector with dim entries. If numeric is true, 0 is inserted
 * at every position. Otherwise a new fraction objects are placed there.
 *
 * @param  {Number} dim
 * The dimension of the vector.
 * 
 * @param  {Boolean} numeric
 * If this parameter is true, the integer 0 is inserted at every position in the
 * vector. Otherwise a new fraction object representing 0 is put there. This
 * parameter is optional, the default value is true. Note that computing an
 * exact result can take significantly more time.
 * 
 * @return {Object}
 * A new null vector with the specified dimension.
 */
function createNullVector(dim, numeric) {
    if (typeof numeric === "undefined") {
        numeric = true;
    }

    var res = [];
    if (numeric) {
        for (var i = 0; i < dim; i++) {
            res[i] = 0;
        }
    } else {
        for (var i = 0; i < dim; i++) {
            res[i] = new Fraction(0, 1);
        }
    }
    return new Vector(res, numeric);
}

/**
 * Returns a new unit vector with dim entries. If numeric is true, 0 is inserted
 * at every position other than n, where the value is 1. Otherwise new fraction
 * objects with the appropriate values are insert.
 *
 * @param  {Number} dim
 * The dimension of the vector.
 *
 * @param {Number} n
 * The positon where to place the 1.
 * 
 * @param  {Boolean} numeric
 * If this parameter is true, the integer 0 is inserted at every position in the
 * vector, except at position n. There the value is 1. Otherwise a new fraction
 * object representing 0 or 1 is put there. This parameter is optional, the
 * default value is true. Note that computing an exact result can take
 * significantly more time.
 * 
 * @return {Object}
 * A new vector with the specified dimension.
 */
function createUnitVector(dim, n, numeric) {
    if (typeof numeric === "undefined") {
        numeric = true;
    }

    var res = [];
    if (numeric) {
        for (var i = 0; i < dim; i++) {
            res[i] = (i === n ? 1 : 0);
        }
    } else {
        for (var i = 0; i < dim; i++) {
            res[i] = (i === n ? new Fraction(1, 1) : new Fraction(0, 1));
        }
    }
    return new Vector(res, numeric);
}

/**
 * Multiplies the fraction with the root and returns the result in simplified
 * form as a string.
 *
 * @param  {Number} frac
 * The fraction.
 * 
 * @param  {Number} root
 * The root.
 * 
 * @return {String}
 * the result in simplified form.
 */
function fracTimesRoot(frac, root) {
    root = root.normalize();
    frac = frac.multiply(new Fraction(root.coef, 1)).normalize();

    if (frac.isZero() || root.rad === 0) {
        return "0";
    } else if (frac.num === 1 && frac.den === 1) {
        return new Root(1, root.rad).toString();
    } else if (root.rad === 1) {
        return frac.toString();
    } else {
        return frac + "*" + new Root(1, root.rad); 
    }
}

/**
 * Returns the position of the first nonnull coefficient (from left to right) or
 * -1 if none is found.
 *
 * @param  {Number[]} equation
 * The equation.
 * 
 * @return {Number}
 * the position of the first nonnull coefficient or -1 if none is found.
 */
function findPivot(equation) {
    for (var i = 0, end = equation.length - 1; i < end; i++) {
        if (equation[i] !== 0) {
            return i;
        }
    }
    return -1;
}

/**
 * Returns an array of the positions of the pivot elements in each equation.
 * The pivot element is the first nonnull coefficient of the equation (read from
 * left to right) or -1 if none is found.
 * 
 * @param  {Number[][]} equations
 * The system of equations.
 *
 * @return {Number[]}
 * the array of the positions of the pivot elements.
 */
function findPivots(equations) {
    var pivots = [];
    for (var i = 0; i < equations.length; i++) {
        pivots[i] = findPivot(equations[i]);
    }
    return pivots;
}

/**
 * Returns a boolean value indicating whether the system of equations is fully
 * specified. That means that it is not empty and that every equation has the
 * same number of elements as any other (an element is either a coefficient or
 * the number on the right hand side). Furthermore each equation must have
 * at least two element. Examples of allowed systems of equations:
 * 
 * 1) [ [1, 2, 3], [1, 2, 3] ]
 * 2) [ [2, 4], [0, 7] ]
 * 
 * Not allowed are the following ones:
 * 
 * 1) null                  (system must have at least one equation)
 * 2) []                    (equation must have at least two elements)
 * 3) [ [1] ]               (same reason)
 * 4) [ [1, 2], [1, 2, 3] ] (number of elements differs)
 * 
 * @param  {Number[][]} equations
 * The system of equations.
 *
 * @return {Boolean}
 * if the system of equations is fully specified.
 */
function isFullySpecified(equations) {
    if (equations.length === 0) {
        return false;
    }
    var len = equations[0].length;
    if (len < 2) {
        return false;
    }
    for (var i = 1; i < equations.length; i++) {
        if (equations[i].length != len) {
            return false;
        }
    }
    return true;
}

/**
 * Creates a string representing the solution for one unknown value. If
 * posible exact values are placed in the string. solution is the vector
 * specifying the solution for this unknown. It has to have the dimension
 * number of unknowns plus 1. The first entry is the coefficient of x0, the
 * next one is the coefficient of x1, and so forth. The rightmost value
 * represents just a number.
 *
 * @param  {Object} solution
 * A vector object specifying the solution for one unknown.
 *
 * @return {String}
 * A string representing the solution for one unknown.
 */
function linSysToStringHelperExact(solution) {
    var elements = solution.elements;
    var isEmpty = true;
    var s = "";
    for (var i = 0, end = solution.dim-1; i <= end; i++) {

        // Coefficient is zero
        if (elements[i].isZero() && (i != end || !isEmpty)) {
            continue;
        }

        // Print operator
        if (!isEmpty) {
            if (elements[i].isNegative()) {
                s += " - ";
            } else {
                s += " + ";
            }
        } else if (elements[i].isNegative()) {
            s += " -";
        }

        isEmpty = false;

        // Coefficient
        if (i !== end && abs(elements[i].num / elements[i].den) !== 1) {
            s += elements[i].abs() + "*";
        } else if (i === end) {
            s += elements[i].abs();
        }

        // Name of unknown
        if (i !== end) {
            s += "x" + i; 
        }
    }
    return s + "\n";
}

/**
 * Creates a string representing the solution for one unknown value. Approximate
 * numeric values are placed in the string. solution is the vector specifying
 * the solution for this unknown. It has to have the dimension
 * number of unknowns plus 1. The first entry is the coefficient of x0, the
 * next one is the coefficient of x1, and so forth. The rightmost value
 * represents just a number.
 *
 * @param  {Object} solution
 * A vector object specifying the solution for one unknown.
 *
 * @return {String}
 * A string representing the solution for one unknown.
 */
function linSysToStringHelperNumeric(solution) {
    var elements = solution.elements;
    var isEmpty = true;
    var s = "";
    for (var i = 0, end = elements.length-1; i <= end; i++) {

        // Coefficient is zero
        if (elements[i] === 0 && (i != end || !isEmpty)) {
            continue;
        }

        // Print operator
        if (!isEmpty) {
            if (elements[i] < 0) {
                s += " - ";
            } else {
                s += " + ";
            }
        } else if (elements[i] < 0) {
            s += " -";
        }

        isEmpty = false;

        // Coefficient
        if (i !== end && abs(elements[i]) !== 1) {
            s += abs(elements[i]) + "*";
        } else if (i === end) {
            s += abs(elements[i]);
        }

        // Name of unknown
        if (i !== end) {
            s += "x" + i; 
        }
    }
    return s + "\n";
}

/**
 * Returns a vector specifying the concrete solution for one unknown value.
 *
 * @param  {Number[]} equation
 * The current equation.
 *
 * @param  {Number} pivot
 * The index of the unknown.
 * 
 * @param  {Object[]} solution
 * An array of the solution vectors.
 *
 * @param  {Boolean} numeric
 * If the solution should be evaluated approximately or not. Note that
 * computing exact results can take significantly more time.
 *
 * @return {Object}
 * A vector representing the solution for the current unknown.
 */ 
function solveLinearSystemCurrent(equation, pivot, solution, numeric) {
    var len = equation.length;
    var current = createNullVector(len, numeric);
    if (numeric) {
        current.elements[len-1] = equation[len-1];
        for (var i = pivot+1; i < len-1; i++) {
            var temp = solution[i].multiply(-equation[i]);
            current = current.add(temp);
        }
        return current.multiply(1 / equation[pivot]);
    } else {
        current.elements[len-1] = new Fraction(equation[len-1], 1);
        for (var i = pivot+1; i < len-1; i++) {
            var temp = solution[i].multiply(new Fraction(-equation[i], 1));
            current = current.add(temp);
        }
        return current.multiply(new Fraction(1, equation[pivot]));
    }
}

/**
 * Computes solution vectors for the unknowns in the equation.
 *
 * @param  {[Number[][]]} equations
 * The system of equations. This has to be in the simplified form returned by
 * {@code gaussianElimination}. For more information refer to the documentation
 * of {@code solveLienarSystem}.
 *
 * @param  {Boolean} numeric
 * If the solution should be approximate or exact. Note that computing exact
 * results can take significantly more time.
 *
 * @return {Object[]}
 * An array of vectors representing the solutions for the unknowns in the
 * equation.
 */
function solveLinearSystemHelper(equations, numeric) {
    var numberElements = equations[0].length;
    var numberEquations = equations.length;
    var pivots = findPivots(equations);
    pivots[numberEquations] = numberElements - 1; // sentinel value
    var solution = [];

    for (var i = numberEquations - 1; i >= 0; i--) {

        // All coefficients are zero
        if (pivots[i] === -1) {
            if (equations[i][numberElements - 1] !== 0) {
                return "This system of equations does not have a solution.";
            } else {
                continue; // No restriction of the solution
            }
        }

        // Unrestricted variables
        for (var j = pivots[i+1] - 1; j > pivots[i]; j--) {
            solution[j] = createUnitVector(numberElements, j, numeric);
        }

        // Current solution
        solution[i] = solveLinearSystemCurrent(
            equations[i],
            pivots[i],
            solution,
            numeric
        );
    }

    // Unrestricted variables
    for (var i = pivots[0] - 1; i >= 0; i--) {
        solution[i] = createUnitVector(numberElements, i, numeric);
    }

    return solution;
}

// =============================================================================
// Math functions
// =============================================================================
 
/**
 * Computes the absolute value of the argument. The absolute value is the
 * argument if it is greater than or equal to zero. Otherwise it is the
 * additive inverse.
 *
 * @param  {Number} n
 * The argument.
 *
 * @return {Number}
 * the absolute value of the argument. 
 */
function abs(n) {
    return n >= 0 ? n : -n;
}

/**
 * Computes the prime factorization of n and returns the result as an array.
 * Each array element is another array with two elements. The first element
 * is the prime and the second one is the exponent. n has to be an integer, if
 * it is negative its absolue value is factorized.
 *
 * @param  {Number} n
 * The number to be factorized. This has to be an integer.
 *
 * @return {Number[][]}
 * the prime factorization of the argument.
 */
function factorize(n) {
    n = abs(n);
    var res = [];
    var insertPos = 0;

    for (var i = 2; n > 1; i++) {
        if (!isPrime(i)) {
            continue;
        }
        var exp = 0;
        while (n > 1 && n % i === 0) {
            exp++;
            n /= i;
        }
        if (exp > 0) {
            res[insertPos] = [i, exp];
            insertPos++;
        }
    }

    return res;
}

/**
 * Implements the Gaussian Elimination. equations is a two dimensional array of
 * equations that has to be fully specified. This is not checked by this
 * function. The content of the array gets changed, so use defensive copying if
 * you need it afterwards. Each element in the array represents one equation.
 *
 * It should be noted, that this function does not solve the equation, but that
 * is simplifies it, so that it can be solved easier. Use 
 * {@code solveLinearSystem} to get the concrete solution.
 * 
 * @param  {Number[][]} equations
 * The system of equations.
 *
 * @param {Boolean} numeric
 * If the solution should be approximate or exact. Note that computing exact
 * results can take significantly more time.
 */
function gaussianElimination(equations, numeric) {
    var numberElements = equations[0].length;
    var numberEquations = equations.length;
    for (var i = 0; i < numberEquations; i++) {
        var c1 = equations[i][i];

        // Swap equations
        var j = i+1;
        while (c1 === 0 && j < numberEquations) {
            var temp = equations[i];
            equations[i] = equations[j];
            equations[j] = temp;
            c1 = equations[i][i];
            j++;
        }
        if (c1 === 0) {
            continue; // Variable eliminated in every equation
        }

        for (var j = i+1; j < numberEquations; j++) {
            var c2 = equations[j][i];
            if (c2 === 0) {
                continue; // Variable is already eliminated
            }
            for (var k = i; k < numberElements; k++) {
                equations[j][k] = c2 * equations[i][k] - c1 * equations[j][k];
            }
        }
    }
}

/**
 * Computes the greatest common divisor of a and b. That is the greatest
 * natural number (including zero) that evenly divides both a and b. Both a and
 * b have to be integers. This is not checked by the function, but has to be
 * enforced by the caller to get correct results.
 *
 * @param  {Number} a
 * The first argument. This has to be an integer.
 *
 * @param  {Number} b
 * The second argument. This has to be an integer.
 *
 * @return {Number}
 * the greatest common divisor of a and b.
 */
function gcd(a, b) {
    return b === 0 ? abs(a) : gcd(b, a % b);
}

/**
 * Returns true if and only if n is a prime number. A prime number is a number
 * that is greater than or equal to two and that has exactly two divisors (1 and
 * the number itself). The argument n has to be an integer. The function does
 * not check that, so in order to get correct result, the caller has to test
 * that.
 *
 * @param  {Number}  n
 * The number to be tested.
 *
 * @return {Boolean}
 * if the argument is prime.
 */
function isPrime(n) {
    if (n === 2 || n === 3) {
        return true;
    }
    if (n < 2 || n % 2 === 0 || n % 3 === 0) {
        return false;
    }
    for (var i = 6, end = Math.floor(Math.sqrt(n)); i <= end; i += 6) {
        if (n % (i-1) === 0 || n % (i+1) === 0) {
            return false;
        }
    }
    return true;
}

/**
 * Returns the maximum of a and b. If a is greater than b, the result is a,
 * otherwise it is b.
 *
 * @param  {Number} a
 * the first argument.
 *
 * @param  {Number} b
 * the second argument.
 *
 * @return {Number}
 * the maximum of a and b.
 */
function max(a, b) {
    return a > b ? a : b;
}

/**
 * Returns the maximum of all elements in the array. That is the element that
 * is greater than or equal to any other element in the array.
 *
 * @param  {Number[]} arr
 * the array of arguments.
 *
 * @return {Number}
 * the maximum of all elements in the array.
 */
function maxArr(arr) {
    if (arr.length === 0) {
        throw new Error("Array must have at least one element.");
    }
    var res = arr[0];
    for (var i = 1; i < arr.length; i++) {
        if (arr[i] > res) {
            res = arr[i];
        }
    }
    return res;
}

/**
 * Returns the minimum of a and b. If a is greater than b, the result is b,
 * otherwise it is a.
 *
 * @param  {Number} a
 * the first argument.
 *
 * @param  {Number} b
 * the second argument.
 *
 * @return {Number}
 * the minimum of a and b.
 */
function min(a, b) {
    return a > b ? b : a;
}

/**
 * Returns the minimum of all elements in the array. That is the element that
 * is less than or equal to any other element in the array.
 *
 * @param  {Number[]} arr
 * the array of arguments.
 *
 * @return {Number}
 * the minimum of all elements in the array.
 */
function minArr(arr) {
    if (arr.length === 0) {
        throw new Error("Array must have at least one element.");
    }
    var res = arr[0];
    for (var i = 1; i < arr.length; i++) {
        if (arr[i] < res) {
            res = arr[i];
        }
    }
    return res;
}

/**
 * Returns a^b. b has to be a natural number or zero, a can be an arbitrary
 * number. This is not enforced by the function, but has to be ensured by the
 * caller.
 *
 * @param  {Number} a
 * Corresponds to the a in a^b.
 *
 * @param  {Number} b
 * Corresponds to the b in a^b. This has to be a natural number or zero.
 *
 * @return {Number}
 * the result of a^b.
 */
function power(a, b) {
    if (b === 0) {
        return 1;
    } else if (b % 2 === 0) {
        var temp = power(a, b / 2);
        return temp * temp;
    } else {
        return a * power(a, b-1);
    }
}

/**
 * Returns the sign of the argument. If it is negative, -1 is the result, if
 * it is zero, the method returns 0, and if it is greater than zero, 1 is
 * returned.
 *
 * @param  {Number} n
 * The argument.
 *
 * @return {Number}
 * the sign of the argument.
 */
function sign(n) {
    if (n < 0) {
        return -1;
    } else if (n === 0) {
        return 0;
    } else {
        return 1;
    }
}

// =============================================================================
// Fraction
// =============================================================================

/**
 * Constructs a new fraction object. num is the numerator of the fraction and
 * den is its denominator. The denominator must not be zero or an error will
 * be thrown.
 * Numerator and denominator are converted to integers. If they have to many
 * decimal places an overflow may occur. No error is thrown if this happens.
 *
 * @param {Number} num
 * The numerator of the fraction.
 * 
 * @param {Number} den
 * The numerator of the fraction. This must be nonzero.
 */
function Fraction(num, den) {
    if (den === 0) {
        throw new Error("Division by zero occured.");
    }

    /**
     * Returns the string s with n zeros appended.
     *
     * @param  {String} s
     * The string.
     * 
     * @param  {Number} n
     * The number of zeros to append.
     * 
     * @return {String}
     * the string with n zeros appended.
     */
    var appendZeros = function(s, n) {
        for (var i = 0; i < n; i++) {
            s += "0";
        }
        return s;
    };

    /**
     * Returns the string a without trailing zeros.
     *
     * @param  {String} s
     * The string.
     * 
     * @return {String}
     * The string without trailing zeros.
     */
    var truncateZeros = function(s) {
        for (var i = s.length - 1; i >= 0; i--) {
            if (s.charAt(i) !== "0") {
                return s.slice(0, i + 1);
            }
        }
        return s;
    };

    var numParts = (num + ".").split(".", 2);
    numParts[1] = truncateZeros(numParts[1]);
    
    var denParts = (den + ".").split(".", 2);
    denParts[1] = truncateZeros(denParts[1]);

    var maxLen = max(numParts[1].length, denParts[1].length);
    numParts[1] = appendZeros(numParts[1], maxLen - numParts[1].length);
    denParts[1] = appendZeros(denParts[1], maxLen - denParts[1].length);

    /**
     * The numerator of the fraction.
     */
    this.num = parseInt(numParts[0] + numParts[1], 10);

    /**
     * The denominator of the fraction.
     */
    this.den = parseInt(denParts[0] + denParts[1], 10);
}

 /**
 * Computes the absolute value of the fraction. The absolute value is the
 * fraction if it is greater than or equal to zero. Otherwise it is the
 * additive inverse.
 *
 * @return {Object}
 * the absolute value of the fraction as a new fraction object. 
 */
Fraction.prototype.abs = function() {
    return new Fraction(abs(this.num), abs(this.den));
};

/**
 * Adds frac to this fraction and returns the result as a new fraction
 * object.
 *
 * @param  {Object} frac
 * The summand. This has to have the properties num and den.
 *
 * @return {Object}
 * the sum of this object and frac.
 */
Fraction.prototype.add = function(frac) {
    var numNew = this.num * frac.den + frac.num * this.den;
    var denNew = this.den * frac.den;
    return new Fraction(numNew, denNew);
};

/**
 * Divides this fraction by frac and returns the result as a new fraction
 * object.
 *
 * @param  {Object} frac
 * The divisor. This has to have the properties num and den.
 * 
 * @return {Object}
 * the quotient of this object and frac.
 */
Fraction.prototype.divide = function(frac) {
    var numNew = this.num * frac.den;
    var denNew = this.den * frac.num;
    return new Fraction(numNew, denNew);
};

/**
 * Returns true if and only if the fraction is less than zero.
 * 
 * @return {Boolean}
 * if this fraction is less than zero.
 */
Fraction.prototype.isNegative = function() {
    return sign(this.num) !== sign(this.den) && !this.isZero();
};

/**
 * Computes the multiplicative inverse of this fraction and returns the result
 * as a new fraction object.
 *
 * @return {Object}
 * the multiplicative inverse of this fraction.
 */
Fraction.prototype.invert = function() {
    if (this.num === 0) {
        throw new Error("Division by zero occured.");
    }
    return new Fraction(this.den, this.num);
};

/**
 * Returns true if and only if the fraction is greater than zero.
 * 
 * @return {Boolean}
 * if this fraction is greater than zero.
 */
Fraction.prototype.isPositive = function() {
    return sign(this.num) === sign(this.den) && !this.isZero();
};

/**
 * Returns true if and only if the fraction equals zero.
 * 
 * @return {Boolean}
 * if this fraction is equal to zero.
 */
Fraction.prototype.isZero = function() {
    return this.num === 0;
};

/**
 * Multiplies this fraction with frac and returns the result as a new fraction
 * object.
 *
 * @param  {Object} frac
 * The factor. This has to have the properties num and den.
 * 
 * @return {Object}
 * the product of this object and frac.
 */
Fraction.prototype.multiply = function(frac) {
    var numNew = this.num * frac.num;
    var denNew = this.den * frac.den;
    return new Fraction(numNew, denNew);
};

/**
 * Computes the additive inverse of this fraction and returns the result as a
 * new fraction object.
 *
 * @return {Object}
 * the additive inverse of this fraction.
 */
Fraction.prototype.negate = function() {
    return new Fraction(-this.num, this.den);
};

/**
 * Simplifies the fraction and normalizes the sign. That means that any
 * superfluos sign in numerator and denominator is removed. If the fraction
 * is negative, the negation is done in the numerator. This method returns
 * a new fraction object. Examples:
 * 
 *  1/-2  --> -1 /2
 * -1/ 2  --> -1 /2
 * -1/-2  -->  1 /2
 *  1/ 2  -->  1 /2
 *  2/-1  -->  -2/1
 *  8/10  -->  4 /5
 *
 * @return {Object}
 * the normalized fraction.
 */
Fraction.prototype.normalize = function() {
    var gcdRes = gcd(this.num, this.den);
    var numNew = this.num / gcdRes;
    var denNew = this.den / gcdRes;
    if (sign(numNew) === sign(denNew)) {
        return new Fraction(abs(numNew), abs(denNew));
    } else {
        return new Fraction(-abs(numNew), abs(denNew));
    }
};

/**
 * Returns the sign of the fraction. If it is negative, -1 is the result, if
 * it is zero, the method returns 0, and if it is greater than zero, 1 is
 * returned.
 *
 * @return {Number}
 * the sign of the fraction.
 */
Fraction.prototype.sign = function() {
    if (isNegative() && !isZero()) {
        return -1;
    } else if (isZero()) {
        return 0;
    } else {
        return 1;
    }
};

/**
 * Subtracts frac from this fraction and returns the result as a new
 * fraction object.
 *
 * @param  {Object} frac
 * The subtrahend. This has to have the properties num and den.
 * 
 * @return {Object}
 * the difference of this object and frac.
 */
Fraction.prototype.subtract = function(frac) {
    var numNew = this.num * frac.den - frac.num * this.den;
    var denNew = this.den * frac.den;
    return new Fraction(numNew, denNew);
};

/**
 * Returns a string representation of the fraction. It is simplified as much
 * as possible. Examples:
 * 
 *  1/-2  --> -1/2
 * -1/ 2  --> -1/2
 * -1/-2  -->  1/2
 *  1/ 2  -->  1/2
 *  0/ 1  -->  0
 *  2/ 1  -->  2
 *  2/-1  -->  -2
 *  8/10  -->  4/5
 *
 * @return {String}
 * a string representation of the simplified fraction.
 */
Fraction.prototype.toString = function() {
    if (this.num % this.den === 0) {
        return this.num / this.den;
    } else {
        var normalized = this.normalize();
        return normalized.num + "/" + normalized.den;
    }
};

// =============================================================================
// Root
// =============================================================================

/**
 * Constructs a new root object. The radicand can be less than zero. In this
 * case it is treated as a complex number.
 *
 * @param {Number} coef
 * The coefficient of the root ( coef*sqrt(rad) ).
 * 
 * @param {[type]} rad
 * The radicand of the root ( coef*sqrt(rad) ).
 */
function Root(coef, rad) {

    /**
     * The coefficient of the root.
     */
    this.coef = coef;

    /**
     * The radicand of the root.
     */
    this.rad = rad;
}

/**
 * Simplifies the root as much as possible. Examples:
 *
 * Root(1, 2) -> Root(1, 2)
 * Root(2, 4) -> Root(4, 1)
 *
 * @return {Object}
 * the simplified root.
 */
Root.prototype.normalize = function() {
    var primeFactors = factorize(this.rad);
    var coefNew = this.coef;
    var radNew = this.rad;
    for (var i = 0; i < primeFactors.length; i++) {
        var prime = primeFactors[i][0];
        var exp = Math.floor(primeFactors[i][1] / 2);
        var powerRes = power(prime, exp);
        coefNew *= powerRes;
        radNew /= powerRes * powerRes;
    }
    return new Root(coefNew, radNew);
};

/**
 * Returns a string representation of this root. It is simplified as much as
 * possible. Examples:
 *
 * Root(1,  2) -> sqrt(2)
 * Root(1, -2) -> sqrt(2)i
 * Root(2,  4) -> 4
 * Root(2,  0) -> 0
 * Root(0,  2) -> 0
 * Root(3,  1) -> 3
 *
 * @return {String}
 * the string representation of the simplified root.
 */
Root.prototype.toString = function() {
    var normalized = this.normalize();
    if (normalized.coef === 0 || normalized.rad === 0) {
        return "0";
    } else if (normalized.coef === 1 && normalized.rad === 1) {
        return "1";
    } else if (normalized.coef === -1 && normalized.rad === 1) {
        return "-1";
    }

    var res = "";

    // Coefficient
    if (normalized.coef === -1) {
        res += "-";
    } else if (normalized.coef !== 1) {
        res += normalized.coef;
    }

    // Radicand
    if (normalized.rad === -1) {
        return res + "i";
    } else if (normalized.rad < 0) {
       return "sqrt(" + (-normalized.rad) + ")i";  
    } else if (normalized.rad !== 1) {
        return "sqrt(" + normalized.rad + ")";  
    }

    return res;
};

// =============================================================================
// Vector
// =============================================================================

/**
 * Constructs a new vector object.
 *
 * @param {Number[] or Object[]} elements
 * Represents the elements of the vector. These can be either numbers or
 * fractions; numeric has to be set accordingly.
 * 
 * @param {Boolean} numeric
 * If operations should be done approximately or exactly. Note that computing
 * exact result can take significantly more time.
 */
function Vector(elements, numeric) {
    if (typeof numeric == "undefined") {
        numeric = true;
    }

    this.elements = elements;
    this.dim = elements.length;
    this.numeric = numeric;
}

/**
 * Adds vec to this vector and returns the result as a new vector. The
 * dimensions and types (numeric or not) of the vectors have to match, otherwise
 * an exception is thrown.
 *
 * @param  {Object} vec
 * The summand as a vector object. This has to have the same dimension and type
 * as this vector.
 *
 * @return {Object}
 * the result as a new vector object.
 */
Vector.prototype.add = function(vec) {
    if (this.dim != vec.dim) {
        throw new Error("Dimensions of vectors do not match.");
    }
    if (this.numeric !== vec.numeric) {
        throw new Error("Types of vectors do not match.");
    }

    var elementsNew = [];
    if (this.numeric) {
        for (var i = 0; i < this.dim; i++) {
            elementsNew[i] = this.elements[i] + vec.elements[i];
        }
    } else {
        for (var i = 0; i < this.dim; i++) {
            elementsNew[i] = this.elements[i].add(vec.elements[i]);
        }
    }
    return new Vector(elementsNew, this.numeric);
};

/**
 * Multiplies vec with this vector and returns the result as a new vector. The
 * dimensions and types (numeric or not) of the vectors have to match, otherwise
 * an exception is thrown.
 *
 * @param  {Object} vec
 * The factor as a vector object. This has to have the same dimension and type
 * as this vector.
 *
 * @return {Object}
 * the result as a new vector object.
 */
Vector.prototype.multiply = function(n) {
    if ((this.numeric && (n instanceof Fraction)) ||
            (!this.numeric && typeof n === "Number")) {
        throw new Error("Types of vector and scalar do not match.");
    }

    var elementsNew = [];
    if (this.numeric) {
        for (var i = 0; i < this.dim; i++) {
            elementsNew[i] = this.elements[i] * n;
        }
    } else {
        for (var i = 0; i < this.dim; i++) {
            elementsNew[i] = this.elements[i].multiply(n);
        }
    }
    return new Vector(elementsNew, this.numeric);
};