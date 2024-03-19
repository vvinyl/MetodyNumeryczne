#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
using namespace std;
typedef double** matrix;
void pamiecNew(matrix& b, const int& m_, const int& n_)
{
    b = new double* [m_];
    for (int i = 0; i < m_; ++i)
        b[i] = new double[n_];
}




void pamiecDel(matrix& b, const int& m_, const int& n_)
{
    for (int i = 0; i < m_; ++i)
        delete[]b[i];

    delete[]b;
}
int wypelnij(matrix& b, const int& m_, const int& n_, const int& w)
{
    for (int i = 0; i < m_; ++i)
        for (int j = 0; j < n_; ++j)
            b[i][j] = w;
    return 0;
}
int wyswietl(const matrix& c, const int& m_, const int& n_)
{
    for (int i = 0; i < m_; ++i)
    {
        for (int j = 0; j < n_; ++j)
            cout << "\t" << c[i][j];
        cout << endl;
    }
    return 0;
}
/*int wypelnij_los(matrix& b, const int& m_, const int& n_, const int& minZ, const int& maxZ)
{
    int zakres = maxZ - minZ + 1;
    srand(time(NULL));
    for (int i = 0; i < m_; ++i)
        for (int j = 0; j < n_; ++j)
            b[i][j] = rand() % zakres + minZ;
    return 0;
}*/
int zapiszDoPliku(string nazwaPliku, const matrix& c, const int& m_, const int& n_)
{
    ofstream plik;
    plik.open(nazwaPliku.c_str());
    for (int i = 0; i < m_; ++i)
    {
        for (int j = 0; j < n_; ++j)
            plik << "\t" << c[i][j];
        plik << endl;
    }
    plik.close();
    return 0;
}
/*int wypelnij_los_sym(matrix& b, const int& m_, const int& n_, const int& startZ, const int& endZ)
{
    int zakres = endZ - startZ + 1;
    srand(time(NULL));
    for (int i = 0; i < m_; ++i)
        for (int j = i; j < n_; ++j)
            b[i][j] = b[j][i] = rand() % zakres + startZ;
    return 0;
}*/
bool czySymetryczna(const matrix& c, const int& m_, const int& n_)
{
    for (int i = 0; i < m_; ++i)
    {
        for (int j = i + 1; j < n_; ++j)
            if (c[i][j] != c[j][i])
                return false;
    }
    return true;

}

int czytajZpliku(string nazwaPliku, matrix& b, const int& m_, const int& n_)
{
    ifstream plikA;
    plikA.open(nazwaPliku.c_str());
    for (int i = 0; i < m_; ++i)
        for (int j = 0; j < n_; ++j)
            plikA >> b[i][j];
    plikA.close();
    return 0;
}

int dodawanie(matrix& a, matrix& b, const int& m_, const int& n_, matrix& s) {
    for (int i = 0; i < m_; ++i)
    {
        for (int j = 0; j < n_; ++j) s[i][j] = a[i][j] + b[i][j];
    }
    return 0;
}
int odejmowanieOdSumy(matrix& s, matrix& c, const int& m_, const int& n_, matrix& so) {
    for (int i = 0; i < m_; ++i)
    {
        for (int j = 0; j < n_; ++j) so[i][j] = s[i][j] - c[i][j];
    }
    return 0;
}
int transponowanie(matrix& b, matrix& bt, const int& m_, const int& n_) {
    for (int i = 0; i < n_; i++) {
        for (int j = 0; j < m_; j++) {
            bt[j][i] = b[i][j];
        }
    }
    return 0;
}

int mnozenie(matrix& a, matrix& bt, matrix& mn, const int& m_, const int& n_) {
    for (int i = 0; i < m_; i++) {
        for (int j = 0; j < n_; j++) {
            mn[i][j] = 0;
            for (int k = 0; k < n_ - 1; k++) {
                mn[i][j] = mn[i][j] + (a[i][k] * bt[k][j]);
            }
        }
    }
    return 0;
}
int mnozenie2(matrix& e, matrix& f, matrix& mn2, const int& m_, const int& n_) {
    for (int i = 0; i < m_; i++) {
        for (int j = 0; j < n_; j++) {
            mn2[i][j] = 0;
            for (int k = 0; k < n_; k++) {
                mn2[i][j] = mn2[i][j] + (e[i][k] * f[k][j]);
            }
        }
    }
    return 0;
}

int mnozenie3(matrix& e, matrix& f, matrix& mn3, const int& m_, const int& n_) {
    for (int i = 0; i < m_; i++) {
        for (int j = 0; j < n_; j++) {
            mn3[i][j] = 0;
            for (int k = 0; k < m_; k++) {
                mn3[i][j] = mn3[i][j] + (e[i][k] * f[k][j]);
            }
        }
    }
    return 0;
}

double wyznacznik(const matrix& e, const int& m_, const int& n_)
{
    double det = e[0][0];
    for (int s = 0; s < n_ - 1; s++) {
        for (int i = s + 1; i < n_; i++) {
            for (int j = s + 1; j < n_; j++) {
                e[i][j] = e[i][j] - (e[i][s] / e[s][s]) * e[s][j];
            }
        }
        det = det * e[s + 1][s + 1];
    }
    return det;
}
int odwracanie(const matrix& mWe, const int& m_, const int& n_, const matrix& mWy) {
    matrix A;
    pamiecNew(A, m_, 2 * n_);//tworzymy dodatkowa macierz ktora ma dwa razy wiecej kolumn

    for (int i = 0; i < m_; i++) {
        for (int j = 0; j < n_; j++) {
            A[i][j] = mWe[i][j];
        }
    }

    for (int i = 0; i < m_; i++) {
        for (int j = 0; j < n_; j++) {
            A[i][j + n_] = 0;
        }
    }

    for (int i = 0; i < m_; i++) {
        A[i][i + n_] = 1;
    }

    wyswietl(A, m_, 2 * n_);
    //
    cout << endl << endl;

    for (int i = 0; i < m_; i++) {
        for (int j = 0; j < n_; j++) {
            mWy[i][j] = A[i][j + n_];
        }
    }

    for (int s = 0; s < n_; s++) {
        double c;
        c = A[s][s];
        A[s][s] = A[s][s] - 1;
        for (int j = s + 1; j < 2 * n_; j++) {
            double d;
            d = A[s][j] / c;
            for (int i = 0; i < n_; i++) {
                A[i][j] = A[i][j] - (d * A[i][s]);
            }
        }
    }
    /*for (int s = 0; s<n_; s++) {
        double c = A[s][s];
        A[s][s]=A[s][s]-1;
        for(int j=s+1;j<2*n_;j++){
            double d = A[s][j]/c;
            for(int i=0; i<n_;i++){
                A[i][j]=A[i][j]-(d*A[i][s]);
            }
        }
    }*/

    cout << "Odwrocona macierz:" << endl;

    for (int i = 0; i < n_; i++) {
        for (int j = n_; j < 2 * n_; j++) {
            if (j == n_ + 3)cout << A[i][j] << endl;
            else cout << A[i][j] << "\t";
        }
    }

    return 0;
}
int main()
{
    cout << "<---------- macierz A ---------->" << endl;
    int mA = 4;
    int nA = 3;
    matrix A;
    pamiecNew(A, mA, nA);
    czytajZpliku("macierzA.txt", A, mA, nA);
    wyswietl(A, mA, nA);

    cout << endl;

    cout << "<---------- macierz B ---------->" << endl;
    int mB = 4;
    int nB = 3;
    matrix B;
    pamiecNew(B, mB, nB);
    czytajZpliku("macierzB.txt", B, mB, nB);
    wyswietl(B, mB, nB);

    cout << endl;

    cout << "<---------- macierz C ---------->" << endl;
    int mC = 4;
    int nC = 3;
    matrix C;
    pamiecNew(C, mC, nC);
    czytajZpliku("macierzC.txt", C, mC, nC);
    wyswietl(C, mC, nC);

    cout << endl;

    cout << "<---------- macierz D ---------->" << endl;
    int mD = 4;
    int nD = 4;
    matrix D;
    pamiecNew(D, mD, nD);
    czytajZpliku("macierzD.txt", D, mD, nD);
    wyswietl(D, mD, nD);

    cout << endl;

    cout << "<---------- macierz E ---------->" << endl;
    int mE = 4;
    int nE = 4;
    matrix E;
    pamiecNew(E, mE, nE);
    czytajZpliku("macierzE.txt", E, mE, nE);
    wyswietl(E, mE, nE);

    cout << endl;

    cout << "<---------- macierz F ---------->" << endl;
    int mF = 4;
    int nF = 1;
    matrix F;
    pamiecNew(F, mF, nF);
    czytajZpliku("macierzF.txt", F, mF, nF);
    wyswietl(F, mF, nF);

    cout << "<---------- DODAWANIE ---------->" << endl;
    int mS = 4;
    int nS = 3;
    matrix S;
    pamiecNew(S, mS, nS);

    dodawanie(A, B, mS, nS, S);
    wyswietl(S, mS, nS);

    cout << "<---------- ODEJMOWANIE OD SUMY ---------->" << endl;
    int mOS = 4;
    int nOS = 3;
    matrix OS;
    pamiecNew(OS, mOS, nOS);
    odejmowanieOdSumy(S, C, mOS, nOS, OS);
    wyswietl(OS, mOS, nOS);
    cout << endl;

    cout << "<---------- TRANSPONOWANIE B ---------->" << endl;
    int mBT = nB;
    int nBT = mB;
    matrix BT;
    pamiecNew(BT, mBT, nBT);
    transponowanie(B, BT, mBT, nBT);
    wyswietl(BT, mBT, nBT);

    cout << "<---------- MNOZENIE A x BT ---------->" << endl;
    int mMN = mA;
    int nMN = nBT;
    matrix MN;
    pamiecNew(MN, mMN, nMN);
    mnozenie(A, BT, MN, mMN, nMN);
    wyswietl(MN, mMN, nMN);
    cout << endl;

    cout << "<---------- MNOZENIE E x D ---------->" << endl;
    int mMN2 = mE;
    int nMN2 = nD;
    matrix MN2;
    pamiecNew(MN2, mMN2, nMN2);
    mnozenie2(E, D, MN2, mMN2, nMN2);
    wyswietl(MN2, mMN2, nMN2);
    cout << endl;

    cout << "<---------- MNOZENIE E x F ---------->" << endl;
    int mMN3 = mE;
    int nMN3 = nF;
    matrix MN3;
    pamiecNew(MN3, mMN3, nMN3);
    mnozenie3(E, F, MN3, mMN3, nMN3);
    wyswietl(MN3, mMN3, nMN3);
    cout << endl;

    cout << "<---------- odwracanie macierzy ---------->" << endl;
    int mMO = mE;
    int nMO = nF;
    matrix MO;
    pamiecNew(MO, mE, nE);
    odwracanie(E, mE, nE, MO);

    cout << "<---------- WYZNACZNIK E---------->" << endl;
    cout << "Wyznacznik macierzy wynosi: " << wyznacznik(E, mE, nE) << endl << endl;

    cout << "<---------- WYZNACZNIK D---------->" << endl;
    cout << "Wyznacznik macierzy wynosi: " << wyznacznik(D, mD, nD) << endl << endl;

    return 0;
}
