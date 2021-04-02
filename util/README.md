# Utility Functions
The language-specific files in this directory contain utility methods to solve Kepler's equation.  The function is called mToE.  The function accepts mean anomaly and eccentricity as function parameters.
Copy just this code to any projects needing the method.  For convenience, an accuracy and timing method is included in each "main" function.  But you only need the mToE implementation in your code.
The timing will be slower than in the original test script because each computation of the eccentric anomaly is a separate function calls.  Most libraries needing this function are implemented to call with one value of mean anomaly and eccentricity at a time.
