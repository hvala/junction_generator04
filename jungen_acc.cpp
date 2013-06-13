
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

using namespace std;

// Take the sum over integers from 1 to n
int sigmaNum(int n)
{
    int sum = 0;
    for ( int i = 1 ; i <= n ; i++ )
    {
        sum += i;
    }
    return sum;
}

// calculate a factorial
int factorial (int n)
{
    int fact = n;
    for( int m = n - 1 ; m > 1 ; m-- )
    {
        fact = fact * m;
    }
    return fact;
}

// partial factorial -- calculate a product over n to k
int part_factorial (int n, int k)
{
    int p_fact = n;
    for ( int m = n - 1 ; m > k ; m-- )
    {
        p_fact = p_fact * m;
    }
    return p_fact;
}

// calculate n choose k
int n_choose_k (int n, int k)
{
    int nck = 0;
    nck = part_factorial(n, n-k) / factorial(k);
    return nck;
}

// express an integer as a string from a different base system -- base < 10
vector<int> change_base(int x, int base)
{
    vector<int> newNumber;

    // find the largest base power to go into the x
    int i = 0;

    do
    {
        ++i;
    } while ( pow(base, i) <= x );
    --i;

    // count the number of times the next greatest base power fits into the number
    while ( i >= 0 )
    {
        int digit = 0;
        do
        {
            digit++;
        } while ( pow(base, i) * digit <= x );
        newNumber.push_back(--digit);

        x = x - ( pow(base, i) * digit );

        --i;
    }

    return newNumber;
}

// take a vector of digits representing a number in a base system 'b', and convert it to a decimal number
int convert_to_decimal (vector<int> number, int base)
{
    vector<int>::iterator iterD;
    int decimal = 0;
    for(iterD = number.begin(), int i = 0 ; iterD < number.end() ; ++iterD, ++i)
    {
        decimal = decimal + (*iterD) * pow(base, i);
    }
    return decimal;
}

// convert a vector of ints to a string
string vec_ints_to_string(vector<int> i)
{
    stringstream str;
    vector<int>::iterator iterI;
    for( iterI = i.begin() ; iterI < i.end() ; iterI++ )
    {
        str << (*iterI);
    }

    string key = str.str();
    return key;
}