// Matrices.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <algorithm>

unsigned debug = 0;

///TODO Fix matrix::echelon on 505.  Test matrix::reduce and matrix::operator~  

using namespace std;

//Strict, weak ordering (yes, there should be a comma between "strict" and "weak" because they are both adjectives)
template<class A>
bool lexrefcomp(vector<A*> x, vector<A*> y) {
    unsigned lenx = x.size(), leny = y.size();
    if (lenx < leny) return true;
    if (lenx > leny) return false;
    
    //In the lenx == 0 case, this loop is skipped correctly.
    for (unsigned i = 0; i < lenx; i++) {
        if (*x[i] < *y[i]) return true;
        if (*x[i] == *y[i]) continue;
        else return false;

    }
    return false;
}

template<class A>
A gcd(A a, A b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}

template<class num>
class mector:public vector<num> {
public:
    mector() {};
    
    //Copy constructor from vectors
    mector(const vector<num>& V) {
        for (auto x : V) {
            (*this).push_back(x);
        }
    };

    //Copy constructor from self
    mector(const mector<num>& M) {
        *this = M;
    }
    
    //Avoids out of range errors.  Wraps indecies around modulo the length of the vector.  Default return is 0 when the dimension is 0.
    num & operator[](const int& index) {
        mector<num> temp(*this);
        unsigned dim = temp.size();

        //Checks for an empty mector
        try {
            if (dim == 0) throw 0;
        }
        catch (int) {
            cout << endl << "Empty mector exception:  You attempted to access a value in an empty mector." << endl << endl;
        }
        
        //Never gives an out of range error and wraps around the mector.  Unlike Python, you can choose any integer and the wrapping will continue.
        unsigned safe_index;
        /*
            This is clunky because % does not work as it would in math.
            The final reduction by dim in the false case is necessary because we could get 0 == (-index % dim) for negative multiples of dim.
            After subtracting 0 from dim, we get dim as our index in the mector.  This is one greater than the largest index.
        */
        safe_index = (index >= 0) ? index : (dim - (-index % dim));
        safe_index %= dim;  //Final reduction by dim
                
        typename mector<num>::iterator it = temp.begin();

        it += safe_index;

        return this->at(safe_index);
    }

    //Adds mectors componentwise.  Prints error and returns an empty mector if the dimensions don't match.
    mector<num> operator+(const mector<num>& a) {
        mector<num> lhs(*this), rhs(a), ret;
        unsigned dim = lhs.size();

        try {
            if (lhs.size() != rhs.size()) throw 1;
        }
        catch (int) {
            cout << "Dimension Error:  Attempted to add differently sized mectors." << endl;
        }
    
        for (unsigned i = 0; i < dim; i++) {
            ret.push_back(lhs[i] + rhs[i]); //Error here with []
        }

        return ret;
    }

    //Dot product.  Prints error if dimentions differ.
    num operator*(const mector<num>& rhs) {
        mector<num> lhs(*this);
        num ret = 0;
        unsigned dim = lhs.size();

        ///Fix this try-catch.  It doesn't actually return an error.
        try {
            if (lhs.size() != rhs.size()) throw 1;
        }
        catch (int) {
            //Dimension error
            cout << "Dimension Error:  Attempted to dot mectors with a different number of entries." << endl;
        }

        mector<num> right;
        right = rhs;
        for (int i = 0; i < dim; i++) {
            ret += lhs[i] * right[i];
        }

        return ret;
    }

    //Lexicographic
    bool operator<(const mector<num>& rhs) {
        mector<num> S, B;
        S = *this;
        B = rhs;

        unsigned size = S.size() < B.size() ? S.size() : B.size();
        
        for (unsigned i = 0; i < size; i++) {
            if (S[i] == B[i]) continue;
            if (S[i] < B[i]) return true;
            return false;
        }
        return false;
    }  
};

//Derferences entries in matrix::rows and matrix::cols as a nonmember to matrix.
template<class num>
mector<mector<num>> deepstar(mector<mector<num*>> A) {
    
    mector<mector<num>> ret;    
    unsigned len = A.size();
    if (len == 0) return ret;
    
    unsigned lillen = A[0].size();

    ret.resize(A.size());

    for (unsigned i = 0; i < len; i++) {
        ret[i].resize(lillen);
    }
    for (unsigned i = 0; i < len; i++) {
        for (unsigned j = 0; j < lillen; j++) {
            ret[i][j] = *A[i][j] ;
        }
    }

    return ret;
}

//Returns 0 for the 0 mector.  Gives an error for an empty mector.
template<class num>
num gcd(mector<num> m) {
    unsigned size = m.size();
    try { if (size == 0) throw 1; }
    catch (int) { cout << "Dimension Error:  Tried to take the gcd of an empty mector."; }

    int ret = m[0];
    ret = ret < 0 ? -ret : ret;  //Absolute value function to avoid negative gcds.  May cause problems with certain numeric classes.
    for (unsigned i = 1; i < size; i++) {
        ret = gcd(ret, m[i]); //This will evaluate gcd about once per entry.  Taking pairs ends up being about the same, but there might be an more efficient approach.
        ret = ret < 0 ? -ret : ret;
    }

    return ret;
}

template<class num1, class num>//Right scalar multiplication (commutative).
mector<num> operator*(const mector<num>& lhs, const num1& rhs) {
    mector<num> temp;
    for (auto x : lhs) {
        temp.push_back(x * rhs);
    }

    return temp;
}

template<class num1, class num>//Left scalar multiplication included for commutativity.
mector<num> operator*(const num1& lhs, const mector<num>& rhs) {
    return rhs * lhs;
}

template<class num>//Cross product.  Produces a vector orthogonal to both inputs (this vector may be 0).
mector<num> cross(const mector<num>& lhs, const mector<num>& rhs) {
    int dim = lhs.size();
    try {
        if (dim != 3 || rhs.size() != 3) throw 1;
    }
    catch (int) {
        cout << "Dimension Error: Attempted to cross different size mectors or mectors that did not have 3 entries." << endl;
        }

    mector<num> ret, left, right;
    left = lhs;
    right = rhs;
    //Think of the determinant model for calculating a cross product.  Here i is the index in the resultant vector.
    for (int i = 0; i < dim; i++) {
        ret.push_back( left[i + 1] * right[i + 2] - left[i - 1] * right[i - 2] );
    }

    return ret;
}


//--------------------------------------------------Matrix Class----------------------------------------------------------

template<class num>
class matrix {
public:
    //The base structure for the matrix.  With dimension and entries, we can deduce everything about the matrix.
    mector<num> entries;
    pair<unsigned, unsigned> dimension;

    //Generated from entries.  These can be added, subtracted, multiplied, divided, dotted, and crossed. (Cross product only works in 3D)
    mector<mector<num*>> rows, cols;

    //Sets all rows based on this->entries
    void set_rows() {
        rows.clear();
        unsigned entry = 0;
        for (unsigned i = 0; i < dimension.first; i++) {

            mector<num*> pushme;
            rows.push_back(pushme);
            for (unsigned j = 0; j < dimension.second; j++) {

                rows[i].push_back(& entries[entry++]);
            }
        }
    }

    //Sets all cols based on this->entries
    void set_cols() {
        cols.clear();
        unsigned entry = 0;
        mector<num*> pushme;
        for (unsigned i = 0; i < dimension.second; i++) {
            cols.push_back(pushme);
        }

        for (unsigned j = 0; j < dimension.first; j++) {
            for (unsigned i = 0; i < dimension.second; i++) {

                cols[i].push_back(&entries[entry++]);
            }
        }
    }

    //Modulus for modular matrices
    unsigned mod = 0;

    //Default constructor
    matrix() {};

    //Without this, VS gave red underlines saying that no suitable copy constructor existed.  To my knowledge, this is equivalent to the default copy constructor.
    matrix(const matrix& original) : entries(original.entries), dimension(original.dimension), rows(original.rows), cols(original.cols), mod(original.mod) {}

    //Somewhat general constructor.  Makes all entries n.  Makes modulus m.
    matrix(unsigned left, unsigned right, num n, unsigned m) {
        entries = vector<num> (left * right, n);
        
        set_rows();
        set_cols();

        mod = m;
    }

    matrix(unsigned left, unsigned right, num n) {
        matrix(left, right, n, 0);
    }

    //I decided the identity matrix having fewer inputs would be more important.  The char c never does anything.  It just tells the constructor that you want all 0s instead of identity.
    matrix(unsigned left, unsigned right, const char c) {
        matrix(left, right, 0, 0);
    }

    matrix(unsigned left, unsigned right, mector<num> ents, unsigned m) {
        ents.resize(left * right);
        
        entries = ents;
        dimension.first = left;
        dimension.second = right;
        mod = m;

        set_rows();
        set_cols();
    }
    
    matrix(unsigned left, unsigned right, mector<num> ents) {
        matrix(left, right, ents, 0);
    }

    //Identity matrix (modular)
    matrix(unsigned n, unsigned mod) {
        const unsigned dim = n * n;
        vector<num> temp(dim, 0);
        mector<num> ents(temp);

        for (unsigned i = 0; i < dim; i += n+1) {
            ents[i] = 1;
        }

        matrix(n, n, ents, mod);
    }
    
    matrix(unsigned n) {
        matrix(n, 0);
    }

    //It basically works, but it's pretty primitive.  Consider making this prettier.
    void display() {
        for (auto X : this->rows) {
            cout << '|';
            for (auto x : X) {
                cout << *x << '\t';
            }
            cout << '|' << endl;
        }
    }

    //Is the matrix square?  0x0 counts as square.
    bool is_square() {
        return (dimension.first == dimension.second);
    };

    //If the modulus is 0, then the matrix is treated as if it is not modular.
    bool is_modular() { return mod != 0; };

    //UNTESTED //Returns rows with dereferenced entries
    mector<mector<num>> deref_rows() {
        mector<mector<num>> derefrow;
        for (auto X : rows) {
            mector<num> temp;
            for (auto x : X) {
                temp.push_back(*x);
            }
            derefrow.push_back(temp);
        }
        return derefrow;
    }

    //UNTESTED //Returns cols with dereferenced entries
    mector<mector<num>> deref_cols() {
        mector<mector<num>> derefcol;
        for (auto X : cols) {
            mector<num> temp;
            for (auto x : X) {
                temp.push_back(*x);
            }
            derefcol.push_back(temp);
        }
        return derefcol;
    }


    //Returns the determinant in the underlying numerical type.  Note that 0x0, 1x1, and 2x2 matrices are special.
    num det() {
        if (!is_square()) return 0;
        unsigned dim = dimension.first;

        //0x0 case is special
        if (dim == 0) return 0;

        //1x1 case is special
        if (dim == 1) {
            return entries[0];
        }

        //2x2 case is special
        if (dim == 2) {
            return *rows[0][0] * (*rows[1][1]) - (*rows[1][0]) * (*rows[0][1]);
        }

        //Possible problem here if the num class can't be 0.
        num ret = 0;

        for (int j = 0; j < dim; j++) {
            //Possible problem here if num class can't be 1;
            num product = 1;
            num tcudorp = 1;
            for (int i = 0; i < dim; i++) {
                product *= *rows[i][i + j];
                tcudorp *= *rows[-i][i + j];
            }
            //+= avoided in case num class doesn't have it
            ret = ret + product - tcudorp;
        }
        return ret;
    };

    //Checks inveribility
    bool is_invertible() {
        num x = det();
        if (x == 0) return false;

        //In this case, num should always have an overload for gcd.  In a PID the ideal generated by x and thix->mod would have to be a unit to return 1.  This cannot be defined properly across Euclidean Domains because the Euclidean Function (or Norm) need not be unique.
        if (is_modular()) {
            return (gcd<num>(x, mod) == 1);
        }

        return true;
    };
    
    //Adds matrices entrywise
    matrix operator+(const matrix& a) {
        unsigned n = dimension.first, m = dimension.second;
        try { if (n != a.dimension.first || m != a.dimension.second) throw 1; }
        catch (int) { cout << "Dimension Error:  Tried to add differently sized matrices."; }
        
        mector<num> sums, rhs;
        rhs = a.entries;
        sums = entries + rhs;
        if (is_modular()) {
            try { if (mod != a.mod) throw 1; }
            catch (int) { cout << "Mod Error:  Tried to add matrices with different moduli."; }

            for (auto x : sums) {
                x %= mod;
            }
        }

        //mod will be 0 in any matrix that doesn't explicitly have a modulus.
        return matrix(n, m, sums, mod);
    };

    //Multiplied matrices in the usual way.  Moduli are included and detected.
    matrix operator*(const matrix& rhs) {
        matrix A = *this, B = rhs;
        try { if (A.dimension.second != B.dimension.first) throw 1; }
        catch (int) { cout << "Dimension Error:  Tried to multiply matrices with incompatible dimensions."; }

        mector<num> dots;
        mector<mector<num>> derefrow, derefcol;

        derefrow = A.deref_rows();
        derefcol = B.deref_cols();


        typename mector<mector<num>>::iterator r = derefrow.begin(), c = derefcol.begin(), lastr = derefrow.end(), lastc = derefcol.end();

        if (A.is_modular()) {
            unsigned M = A.mod;
            while (r != lastr) {
                c = derefcol.begin();//Resets the iterator so we get all pairs of rows and cols.

                while (c != lastc) {
                    dots.push_back(((*r) * (*c)) % M);
                    c++; //I love it when this shows up.
                }
                r++;
            }
        }
        else {
            while (r != lastr) {
                c = derefcol.begin();//Resets the iterator so we get all pairs of rows and cols.

                while (c != lastc) {
                    dots.push_back((*r) * (*c));
                    c++; //I love it when this shows up.
                }
                r++;
            }
        }

        //A.mod will be 0 whenever A is not modular.
        return matrix(A.dimension.first, B.dimension.second, dots, A.mod);
    };    

    //Sorts rows for row reduction.
    matrix rowsort() {
        matrix<num> temp(*this);

        sort(temp.rows.begin(), temp.rows.end(), lexrefcomp<num>);

        reverse(temp.rows.begin(), temp.rows.end());
        mector<num> update;

        for (auto X : temp.rows) {
            for (auto x : X) {
                update.push_back(*x);
            }
        }

        return temp;
    }
    
    //BUGGED //FIX NEGATIVE NUMBERS //MOD VERSION UNWRITTEN //Simplifies the matrix into echelon form based on and error value.  Numbers will be rounded to zero when they are less than the error.  Requires a division operator that is equivalent to multiplying by an inverse.
    matrix echelon(const num error) {
        if (dimension.first == 0 || dimension.second == 0) return *this;
      
        matrix<num> temp(*this);

        temp = temp.rowsort();

        if (is_modular()) {
            //Echelon reduction changes quite a lot for modular matrices.  Fill this in later.
           
        }
        else {

            for (unsigned i = 0; i < temp.dimension.first; i++) {
                unsigned pivot = 0;

                //Move pivot position to the right to first nonzero entry.
                while (*temp.rows[i][pivot] == 0 && pivot != temp.dimension.second - 1) { pivot++; }
                

                if (pivot == temp.dimension.second - 1) { 
                    if (temp.entries[-1] != 0) { 
                        //These 2 following lines should be redundant.
                        //*temp.rows[-1][-1] = 1;
                        //*temp.cols[-1][-1] = 1;
                        temp.entries[-1] = 1; 
                    }
                    //Since the rows are sorted, all remaining rows must also be all zeros.
                    return temp; 
                }   
                
                //Possible rounding issues here.
                num inverse = 1 / *temp.rows[i][pivot];
                //*temp.rows[i][pivot] = 1;

                for (unsigned j = pivot; j < temp.dimension.first; j++) {
                    *temp.rows[i][j] = *temp.rows[i][j] * inverse;
                }
                
                //Round based on error.
                for (unsigned i = 0; i < temp.entries.size(); i++) {
                    num k = temp.entries[i];
                    if (k > -error && k < error ) temp.entries[i] = 0;
                }

                /*
                for (unsigned k = 0; k < temp.rows[i].size(); k++) { 
                    temp.rows[i][k] = swap1[k];
                }
                */

                for (unsigned j = i + 1; j < dimension.first; j++) {

                    //swap1 depricated
                    mector<num>swap2(deepstar(temp.rows)[j]);
                    
                    //Set up row subtraction
                    mector<num> rowadd(deepstar(temp.rows)[i]);
                    rowadd = *temp.rows[j][pivot] * rowadd;
                    rowadd = -1 * rowadd;
                    
                    swap2 = swap2 + rowadd;

                    //Update rows
                    for (unsigned k = 0; k < dimension.second; k++) {
                        *temp.rows[j][k] = swap2[k];
                    }
                    //Old version when rows weren't references.
                    //temp.rows[j].swap(swap2);                    
                }

                temp = temp.rowsort();
            }
        }        
        
        if (temp.entries[-1] != 0) {
            //These next 2 lines should be redundant.
            //*temp.rows[-1][-1] = 1;
            //*temp.cols[-1][-1] = 1;
            temp.entries[-1] = 1;
        }

        /* 
        //Old version when rows weren't references
        mector<num> update;
        for (auto X : temp.rows) {
            for (auto x : X) {
                update.push_back(x);
            }
        }
                
        temp = matrix<num>(temp.dimension.first, temp.dimension.second, update, temp.mod);
        */
        
        return temp;
    }
    

    //UNTESTED //Returns the Row Reduced form of the matrix.  Does nothing for modular matrices yet.
    matrix reduce(const num error) {
        if (dimension.first == 0 || dimension.second == 0) return *this;
        matrix ret = *this;
        ret = ret.echelon(error);
        //Removes reverences to simplify mector addition.
        mector<mector<num>> temp;
        for (auto X : ret.rows) {
            mector<num> a;
            temp.push_back(a);
            for (auto *x : X) {
                a.push_back(*x);
            }
        }

        if (ret.is_modular()) {
            //Writeme.  Plz be careful.  This might get weird.
        }
        else {
            //No need to check the first row.
            for (unsigned i = 1; i < dimension.first; i++) {
                //No need to check the first column.
                for (unsigned j = 1; j < dimension.second; j++) {
                    if (temp[i][j] == 1) {
                        for (int k = i - 1; k >= 0; k--) {
                            num clearme = temp[k][j];
                            if (clearme != 0) {
                                temp[k] = temp[k] + (-temp[k][j] * temp[i]);
                            }
                        }
                    }
                }
            }
        }

        /*
        //Old version from when cols weren't references
        mector<num> update;
        for (unsigned y = 0; y < dimension.first; y++) {
            for (unsigned x = 0; x < dimension.second; x++) {
                update.push_back(ret.cols[x][y]);
            }
        }

        ret = matrix(dimension.first, dimension.second, update, mod);
        */
        
        mector<num> flattemp;
        for (auto X : temp) {
            for (auto x : X) {
                flattemp.push_back(x);
            }
        }
        return matrix<num>(dimension.first, dimension.second, flattemp, mod);
    };

    //UNTESTED
    void rpush_back(const mector<num*> row) {
        rows.push_back(row);
        mector<num> temp;
        for (auto X : rows) {
            for (auto x : X) {
                temp.push_back(*x);
            }
        }
        *this = matrix(dimension.first, dimension.second, temp, mod);
    }

    //UNTESTED
    void rpop_back() {
        rows.pop_back();
        mector<num> temp;
        for (auto X : rows) {
            for (auto x : X) {
                temp.push_back(*x);
            }
        }
        *this = matrix(dimension.first, dimension.second, temp, mod);
    }

    //UNTESTED
    void rpop_front() {
        reverse(rows.begin(), rows.end());
        rows.pop_back();
        reverse(rows.begin(), rows.end());
        
        mector<num> temp;
        for (auto X : rows) {
            for (auto x : X) {
                temp.push_back(*x);
            }
        }
        *this = matrix(dimension.first, dimension.second, temp, mod);
    }

    //Returns the transpose of a matrix.
    matrix transpose() {
        mector<num> update;
        for (auto X : cols) {
            for (auto x : X) {
                update.push_back(*x);
            }
        }

        return matrix<num>(dimension.second, dimension.first, update, mod);
    }

    //UNTESTED
    void cpush_back(const mector<num*> col) {
        //Writeme
    }

    //UNTESTED
    void cpop_back() {

    }

    //UNTESTED
    void cpop_front() {

    }

    //UNTESTED  //Throws exception for noninvertible matrices and assumes that num class throws an exception for impossible division. //Uses 0.001 as default error.
    matrix operator~() {
        matrix<num> temp; 
        temp = *this;
        
        if (is_modular()) {
            //Writeme
        }

        //Invertibility check
        try {
            if (!is_invertible()) throw 1;
            num x = 1 / (det());
        }
        catch (...) {
            cout << "Attempted to invert an uninvertible matrix" << endl;
        }

        unsigned n = dimension.first;

        const matrix<num> I(n, mod);
        for (auto col : I.cols) {
            temp.cpush_back(col);
        }

        temp = temp.reduce(0.001);  cout << endl << "Default error 0.001 used for row reduction." << endl << endl;

        for (auto row : I.cols) {
            temp.cpop_front();
        }

        return temp;
    };
};

//Multiplies each entry of a matrix by a scalar (commutative).
template<class num>
matrix<num> operator*(const matrix<num>& A, const num& s) {
    return matrix<num>(A.dimension.first, A.dimension.second, A.entries * s, A.mod);
}

//Included for commutativity.
template<class num>
matrix<num> operator*(const num& s, const matrix<num>& A) {
    return A * s;
}

int main()
{
    ///Matrix tests
    
    vector<float> v(9, 2);
    v[0] = 2;
    mector<float> threes(v);  //It's actually twos though.
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
    v[3] = 0;
    v[4] = 0;
    v[5] = 0;
    v[6] = 0;
    v[7] = 0;
    v[8] = 1;

    mector<float> twos(v);  //This mector is not actually twos.
    matrix<float> square(3, 3, twos, 0);

    mector<float> sixentries(v);
    sixentries.resize(6);
    matrix<float> rect(3, 2, sixentries, 0);
    
    matrix<float>square2(square);
    square2.entries[0] = 50;
    square.display();
    cout << endl;
    
    //square.echelon(0.001).display();
    square.echelon(0.001).display();
    cout << endl;


    //Mector tests
    /*
    mector<int> buddy;
    for (int i = 0; i < 5; i++) {
        buddy.push_back(i);
    }
    buddy[0] = 12;
    buddy[1] = -66;
    for (int i = 0; i < 5; i++) {
        cout << buddy[i] << endl;
    }
    */
    /*
    mector<int> test2, test3;
    test3.push_back(1);
    test3.push_back(2);
    test3.push_back(3);
    
    test2 = cross(test, test3);


    cout << test2 * test << endl;
    cout << test2 * test3 << endl;

    for (auto x : test) {
        cout << x << ' ';
    }
    cout << endl;
    */
    /*
    mector<int>::iterator it = test.begin();
    for (int i = 0; i < 6; i++) {        
        it += i;
        cout << i << " : " << *it << endl;
        it = test.begin();
    }
    
    int dim = test.size();
    for (int i = -2*dim; i < 2*dim; i++) {
        cout << i << " : " << ((i % dim) + dim) % dim << " : " << test[i] << endl;
    }
    cout << endl << endl;
    */
    /*
    for (auto x : test ) {
        cout << x << ' ';
    }
    cout << endl;
    */
    
    
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
