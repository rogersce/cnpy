#include"cnpy.h"
#include<complex>
#include<cstdlib>
#include<iostream>
#include<map>
#include<string>

const int Nx = 128;
const int Ny = 64;
const int Nz = 32;

int main()
{
    //set random seed so that result is reproducible (for testing)
    srand(0);
    //create random data
    std::vector<std::complex<double>> data(Nx*Ny*Nz);
    for(int i = 0;i < Nx*Ny*Nz;i++) data[i] = std::complex<double>(rand(),rand());

    //save it to file
    cnpy::npy_save("arr1.npy",&data[0],{Nz,Ny,Nx},"w");

    //load it into a new array
    cnpy::NpyArray arr = cnpy::npy_load("arr1.npy");
    std::complex<double>* loaded_data = arr.data<std::complex<double>>();
    
    //make sure the loaded data matches the saved data
    assert(arr.word_size == sizeof(std::complex<double>));
    assert(arr.shape.size() == 3 && arr.shape[0] == Nz && arr.shape[1] == Ny && arr.shape[2] == Nx);
    for(int i = 0; i < Nx*Ny*Nz;i++) assert(data[i] == loaded_data[i]);

    //append the same data to file
    //npy array on file now has shape (Nz+Nz,Ny,Nx)
    cnpy::npy_save("arr1.npy",&data[0],{Nz,Ny,Nx},"a");

    //now write to an npz file
    //non-array variables are treated as 1D arrays with 1 element
    double myVar1 = 1.2;
    char myVar2 = 'a';
    cnpy::npz_save("out.npz","myVar1",&myVar1,{1},"w"); //"w" overwrites any existing file
    cnpy::npz_save("out.npz","myVar2",&myVar2,{1},"a"); //"a" appends to the file we created above
    cnpy::npz_save("out.npz","arr1",&data[0],{Nz,Ny,Nx},"a"); //"a" appends to the file we created above

    //load a single var from the npz file
    cnpy::NpyArray arr2 = cnpy::npz_load("out.npz","arr1");

    //load the entire npz file
    cnpy::npz_t my_npz = cnpy::npz_load("out.npz");
    
    //check that the loaded myVar1 matches myVar1
    cnpy::NpyArray arr_mv1 = my_npz["myVar1"];
    double* mv1 = arr_mv1.data<double>();
    assert(arr_mv1.shape.size() == 1 && arr_mv1.shape[0] == 1);
    assert(mv1[0] == myVar1);

    //now test compressed npz files
    std::cout << "Testing compressed npz save/load..." << std::endl;
    
    //save compressed data
    cnpy::npz_save_compressed("out_compressed.npz","myVar1",&myVar1,{1},"w");
    cnpy::npz_save_compressed("out_compressed.npz","arr1",&data[0],{Nz,Ny,Nx},"a");
    
    //load and verify the compressed file
    cnpy::npz_t my_npz_compressed = cnpy::npz_load("out_compressed.npz");
    
    //check myVar1
    cnpy::NpyArray arr_mv1_c = my_npz_compressed["myVar1"];
    double* mv1_c = arr_mv1_c.data<double>();
    assert(arr_mv1_c.shape.size() == 1 && arr_mv1_c.shape[0] == 1);
    assert(mv1_c[0] == myVar1);
    
    //check arr1
    cnpy::NpyArray arr1_c = my_npz_compressed["arr1"];
    std::complex<double>* arr1_data_c = arr1_c.data<std::complex<double>>();
    assert(arr1_c.shape.size() == 3 && arr1_c.shape[0] == Nz && arr1_c.shape[1] == Ny && arr1_c.shape[2] == Nx);
    for(int i = 0; i < Nx*Ny*Nz; i++) assert(data[i] == arr1_data_c[i]);
    
    //also test loading a single variable from compressed file
    cnpy::NpyArray arr1_single = cnpy::npz_load("out_compressed.npz", "arr1");
    std::complex<double>* arr1_single_data = arr1_single.data<std::complex<double>>();
    for(int i = 0; i < Nx*Ny*Nz; i++) assert(data[i] == arr1_single_data[i]);
    
    std::cout << "Compressed npz test passed!" << std::endl;
}
