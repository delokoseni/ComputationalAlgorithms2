#include <iostream>
#include <fstream>
#include <Windows.h>
#include <vector>

void MatrixOutput(std::vector<std::vector<double>> Matrix, int Dimension, std::ostream& Stream);
void VectorOutput(std::vector<double>Vector, int Dimension, std::ostream& Stream);
int JacobisMethod(std::vector<std::vector<double>> Matrix, std::vector<double> FreeTermsColumn,
                  std::vector<double> *const X, int Dimension, double Accuracy);
int SeidelsMethod(std::vector<std::vector<double>> Matrix, std::vector<double> FreeTermsColumn,
                   std::vector<double> *const X, int Dimension, double Accuracy);
double SumOfRowWithoutOneElement(std::vector<std::vector<double>> Matrix, std::vector<double> X, int Dimension, int IndexOfRow);
bool CheckDiagonal(std::vector<std::vector<double>> Matrix, int Dimension);
bool CheckAccuracy(std::vector<double> X, std::vector<double> XPrevious, int Dimension, double Accuracy);
void SetStartingSolution(std::vector<double>* const X, int Dimension);

int main()
{
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    std::string DimensionFileName = "Dimension.txt";
    std::ifstream InputDimensionFile(DimensionFileName);
    int Dimension = 0; 
    if (InputDimensionFile.is_open())
    {
        InputDimensionFile >> Dimension;
        InputDimensionFile.close();
    }
    else
        std::cout << "Ошибка! Не удалось открыть файл для ввода размерности матрицы!" << std::endl;

    std::string MatrixFileName = "A0.txt";
    std::ifstream InputMatrixFile(MatrixFileName);
    std::vector<std::vector<double>> Matrix(Dimension, std::vector<double>(Dimension));
    if (InputMatrixFile.is_open())
    {
        for(int i = 0; i < Dimension; i++)
            for (int j = 0; j < Dimension; j++)
                InputMatrixFile >> Matrix[i][j];
        InputMatrixFile.close();
    }
    else
        std::cout << "Ошибка! Не удалось открыть файл для ввода матрицы!" << std::endl;

    std::string FreeTermsColumnFileName = "B0.txt";
    std::ifstream InputFreeTermsColumnFile(FreeTermsColumnFileName);
    std::vector<double> FreeTermsColumn(Dimension);
    if (InputFreeTermsColumnFile.is_open())
    {
        for (int i = 0; i < Dimension; i++)
            InputFreeTermsColumnFile >> FreeTermsColumn[i];
        InputFreeTermsColumnFile.close();
    }
    else
        std::cout << "Ошибка! Не удалось открыть файл для ввода столбца свободных членов!" << std::endl;

    std::string AccuracyFileName = "Accuracy.txt";
    std::ifstream InputAccuracyFile(AccuracyFileName);
    double Accuracy = 0.0;
    if (InputAccuracyFile.is_open())
    {
        InputAccuracyFile >> Accuracy;
        InputAccuracyFile.close();
    }
    else
        std::cout << "Ошибка! Не удалось открыть файл для ввода точности!" << std::endl;

    std::string OutputFileName = "Output.txt";
    std::ofstream OutputFile(OutputFileName);
    if(OutputFile.is_open())
    {
        OutputFile << "Матрица: " << std::endl;
        MatrixOutput(Matrix, Dimension, OutputFile);
        OutputFile << std::endl << "Столбец свободных членов: " << std::endl;
        VectorOutput(FreeTermsColumn, Dimension, OutputFile);
        OutputFile << std::endl << "Точность: " << std::endl;
        OutputFile << Accuracy << std::endl;
        OutputFile << std::endl << "Количество итераций метода Якоби: " << std::endl;
        std::vector<double> X(Dimension);
        OutputFile << JacobisMethod(Matrix, FreeTermsColumn, &X, Dimension, Accuracy) << std::endl;
        OutputFile << std::endl << "Решения СЛАУ метода Якоби: " << std::endl;
        VectorOutput(X, Dimension, OutputFile);
        OutputFile << std::endl << "Количество итераций метода Зейделя: " << std::endl;
        OutputFile << SeidelsMethod(Matrix, FreeTermsColumn, &X, Dimension, Accuracy) << std::endl;
        OutputFile << std::endl << "Решения СЛАУ метода Зейделя: " << std::endl;
        VectorOutput(X, Dimension, OutputFile);
        OutputFile.close();
    }
    else
        std::cout << "Ошибка! Не удалось открыть файл для вывода!" << std::endl;
    return 0;
}

void MatrixOutput(std::vector<std::vector<double>> Matrix, int Dimension, std::ostream& Stream) 
{
    for (int i = 0; i < Dimension; i++)
    {
        for (int j = 0; j < Dimension; j++)
        {
            Stream << Matrix[i][j] << "\t\t";
        }
        Stream << std::endl;
    }
}

void VectorOutput(std::vector<double>Vector, int Dimension, std::ostream& Stream)
{
    for (int i = 0; i < Dimension; i++)
    {
        Stream << Vector[i] << std::endl;
    }
}

double SumOfRowWithoutOneElement(std::vector<std::vector<double>> Matrix, std::vector<double> X, int Dimension, int IndexOfRow)
{
    double Sum = 0;
    for(int i = 0; i < Dimension; i++)
        if(IndexOfRow != i)
            Sum += Matrix[IndexOfRow][i] * X[i];
    return Sum;
}

int JacobisMethod(std::vector<std::vector<double>> Matrix, std::vector<double> FreeTermsColumn,
                              std::vector<double> *const X, int Dimension, double Accuracy)
{
    SetStartingSolution(X, Dimension);
    int k = 0;
    if(CheckDiagonal(Matrix, Dimension))
    {
        for(k = 0; k > -1; k++)
        {
            std::vector<double> XPrevious = *X;
            VectorOutput(XPrevious, Dimension, std::cout);
            for (int i = 0; i < Dimension; i++)
                (*X)[i] = (FreeTermsColumn[i] - SumOfRowWithoutOneElement(Matrix, XPrevious, Dimension, i)) / Matrix[i][i];
            if (CheckAccuracy(*X, XPrevious, Dimension, Accuracy))
                break;
        }
    }
    return k;
}

int SeidelsMethod(std::vector<std::vector<double>> Matrix, std::vector<double> FreeTermsColumn,
                  std::vector<double>* const X, int Dimension, double Accuracy)
{
    SetStartingSolution(X, Dimension);
    int k = 0;
    if (CheckDiagonal(Matrix, Dimension))
    {
        for (k = 0; k > -1; k++)
        {
            std::vector<double> XPrevious = *X;
            for (int i = 0; i < Dimension; i++)
                (*X)[i] = (FreeTermsColumn[i] - SumOfRowWithoutOneElement(Matrix, *X, Dimension, i)) / Matrix[i][i];
            if (CheckAccuracy(*X, XPrevious, Dimension, Accuracy))
                break;
        }
    }
    return k;
}

bool CheckDiagonal(std::vector<std::vector<double>> Matrix, int Dimension)
{
    for (int i = 0; i < Dimension; i++)
        if (Matrix[i][i] == 0)
            return false;
    return true;
}

bool CheckAccuracy(std::vector<double> X, std::vector<double> XPrevious, int Dimension, double Accuracy)
{
    double max = fabs(X[0] - XPrevious[0]);
    for (int i = 1; i < Dimension; i++)
        if (fabs(X[i] - XPrevious[i]) > max)
            max = fabs(X[i] - XPrevious[i]);
    if (max < Accuracy)
        return true;
    else
        return false;
}

void SetStartingSolution(std::vector<double>* const X, int Dimension)
{
    for (int i = 0; i < Dimension; i++)
    {
        (*X)[i] = 0.0;
    }
}