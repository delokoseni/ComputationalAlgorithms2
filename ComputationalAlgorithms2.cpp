#include <iostream>
#include <fstream>
#include <Windows.h>
#include <vector>

void MatrixOutput(std::vector<std::vector<double>> Matrix, int Dimension, std::ostream& Stream);
void VectorOutput(std::vector<double>Vector, int Dimension, std::ostream& Stream);

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

    std::string MatrixFileName = "Matrix.txt";
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

    std::string FreeTermsColumnFileName = "FreeTermsColumn.txt";
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
    double Accurancy = 0.0;
    if (InputAccuracyFile.is_open())
    {
        InputAccuracyFile >> Accurancy;
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
        OutputFile << Accurancy << std::endl;
        OutputFile.close();
    }
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