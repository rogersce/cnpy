//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#include"cnpy.h"
#include<complex>
#include<cstdlib>
#include<algorithm>
#include<cstring>
#include<iomanip>
#include<stdint.h>
#include<stdexcept>
#include <regex>

char cnpy::BigEndianTest() {
    int x = 1;
    return (((char *)&x)[0]) ? '<' : '>';
}

char cnpy::map_type(const std::type_info& t)
{
    if(t == typeid(float) ) return 'f';
    if(t == typeid(double) ) return 'f';
    if(t == typeid(long double) ) return 'f';

    if(t == typeid(int) ) return 'i';
    if(t == typeid(char) ) return 'i';
    if(t == typeid(short) ) return 'i';
    if(t == typeid(long) ) return 'i';
    if(t == typeid(long long) ) return 'i';

    if(t == typeid(unsigned char) ) return 'u';
    if(t == typeid(unsigned short) ) return 'u';
    if(t == typeid(unsigned long) ) return 'u';
    if(t == typeid(unsigned long long) ) return 'u';
    if(t == typeid(unsigned int) ) return 'u';

    if(t == typeid(bool) ) return 'b';

    if(t == typeid(std::complex<float>) ) return 'c';
    if(t == typeid(std::complex<double>) ) return 'c';
    if(t == typeid(std::complex<long double>) ) return 'c';

    else return '?';
}

template<> std::vector<char>& cnpy::operator+=(std::vector<char>& lhs, const std::string rhs) {
    lhs.insert(lhs.end(),rhs.begin(),rhs.end());
    return lhs;
}

template<> std::vector<char>& cnpy::operator+=(std::vector<char>& lhs, const char* rhs) {
    //write in little endian
    size_t len = strlen(rhs);
    lhs.reserve(len);
    for(size_t byte = 0; byte < len; byte++) {
        lhs.push_back(rhs[byte]);
    }
    return lhs;
}

void cnpy::parse_npy_header(std::string &buffer, size_t &word_size,
                            std::vector<size_t> &shape, bool &fortran_order) {
    assert(buffer.substr(0, 6) == "\x93NUMPY");
    uint8_t major_version = *reinterpret_cast<uint8_t *>(buffer.data() + 6);
    uint8_t minor_version = *reinterpret_cast<uint8_t *>(buffer.data() + 7);
    uint16_t header_len =
      (*reinterpret_cast<uint8_t *>(buffer.data() + 8)) |
      ((*reinterpret_cast<uint8_t *>(buffer.data() + 9)) << 8);
    std::string header = buffer.substr(10, header_len);

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order")+16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stoi(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable. 
    //not sure when this applies except for byte array
    loc1 = header.find("descr")+9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);

    //char type = header[loc1+1];
    //assert(type == map_type(T));

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}

void cnpy::parse_npy_header(FILE* fp, size_t& word_size, std::vector<size_t>& shape, bool& fortran_order) {  
    char buffer[256];
    size_t res = fread(buffer,sizeof(char),11,fp);       
    if(res != 11)
        throw std::runtime_error("parse_npy_header: failed fread");
    std::string header = fgets(buffer,256,fp);
    assert(header[header.size()-1] == '\n');

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'fortran_order'");
    loc1 += 16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");
    if (loc1 == std::string::npos || loc2 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: '(' or ')'");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stoi(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable. 
    //not sure when this applies except for byte array
    loc1 = header.find("descr");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'descr'");
    loc1 += 9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);

    //char type = header[loc1+1];
    //assert(type == map_type(T));

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}

void cnpy::parse_zip_footer(FILE* fp, uint16_t& nrecs, size_t& global_header_size, size_t& global_header_offset)
{
    std::vector<char> footer(22);
    fseek(fp,-22,SEEK_END);
    size_t res = fread(&footer[0],sizeof(char),22,fp);
    if(res != 22)
        throw std::runtime_error("parse_zip_footer: failed fread");

    uint16_t disk_no, disk_start, nrecs_on_disk, comment_len;
    disk_no = *(uint16_t*) &footer[4];
    disk_start = *(uint16_t*) &footer[6];
    nrecs_on_disk = *(uint16_t*) &footer[8];
    nrecs = *(uint16_t*) &footer[10];
    global_header_size = *(uint32_t*) &footer[12];
    global_header_offset = *(uint32_t*) &footer[16];
    comment_len = *(uint16_t*) &footer[20];

    assert(disk_no == 0);
    assert(disk_start == 0);
    assert(nrecs_on_disk == nrecs);
    assert(comment_len == 0);
}

cnpy::NpyArray load_the_npy_file(FILE* fp) {
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    cnpy::parse_npy_header(fp,word_size,shape,fortran_order);

    cnpy::NpyArray arr(shape, word_size, fortran_order);
    size_t nread = fread(arr.data<char>(),1,arr.num_bytes(),fp);
    if(nread != arr.num_bytes())
        throw std::runtime_error("load_the_npy_file: failed fread");
    return arr;
}

cnpy::NpyArray load_the_npz_array(libzip::file &file, uint64_t size) {
  std::string buffer = file.read(size);

  std::vector<size_t> shape;
  size_t word_size;
  bool fortran_order;
  cnpy::parse_npy_header(buffer, word_size,
                         shape, fortran_order);

  cnpy::NpyArray array(shape, word_size, fortran_order);

  uint64_t offset = size - array.num_bytes();
  memcpy(array.data<unsigned char>(), buffer.data() + offset,
         array.num_bytes());

  return array;
}

cnpy::npz_t cnpy::npz_load(std::string fname) {
  libzip::archive zip(fname, ZIP_RDONLY);
  cnpy::npz_t arrays;

  for (auto stat : zip) {
    if (!(stat.valid & (ZIP_STAT_NAME | ZIP_STAT_SIZE))) {
      // Skip files without name or size
      continue;
    }
    std::string name(stat.name);
    auto file = zip.open(name);

    // erase the lagging .npy
    std::string varname = name.substr(0, name.size() - 4);
    arrays[varname] = load_the_npz_array(file, stat.size);

  }

  return arrays;
}

cnpy::NpyArray cnpy::npz_load(std::string fname, std::string varname) {
  libzip::archive zip(fname, ZIP_RDONLY);
  std::string varfilename = varname + ".npy";
  libzip::stat stat = zip.stat(varfilename);
  if (!(stat.valid & ZIP_STAT_SIZE)) {
    throw std::runtime_error("npz_load: No size associated with " + varname);
  }

  auto file = zip.open(varfilename);
  return load_the_npz_array(file, stat.size);
}

cnpy::NpyArray cnpy::npy_load(std::string fname) {

    FILE* fp = fopen(fname.c_str(), "rb");

    if(!fp) throw std::runtime_error("npy_load: Unable to open file "+fname);

    NpyArray arr = load_the_npy_file(fp);

    fclose(fp);
    return arr;
}



