# ifndef ADMIN_HPP
# define ADMIN_HPP

# include <map>
# include <cmath>
# include <string>
# include <chrono>
# include <vector>
# include <complex>
# include <fftw3.h>
# include <sstream>
# include <iomanip>
# include <fstream>
# include <iostream>
# include <algorithm>

int nextpow2(int n); 

bool str2bool(std::string s);

std::string format1Decimal(float x);

void import_binary_float(std::string path, float * array, int n);
void export_binary_float(std::string path, float * array, int n);

void import_text_file(std::string path, std::vector<std::string> &elements);
 
std::string catch_parameter(std::string target, std::string file);

std::vector<std::string> split(std::string s, char delimiter);

# endif
