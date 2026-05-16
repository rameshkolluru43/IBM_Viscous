#ifndef SOLVERFUNCTIONS_H
#define SOLVERFUNCTIONS_H

/**
 * @brief Generic initializer function
 * @param a1 First value
 * @param a2 Second value
 * @return Maximum of the two values
 */
template <typename T>
T Initialiser(T a1, T a2)
{
	if (a1 >= a2)
		return a1;
	else
		return a2;
}

#endif // SOLVERFUNCTIONS_H
