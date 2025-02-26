#pragma once

#include <iostream>
#include <string>


namespace hydro
{

void initialize(int& argc, char**& argv);
void abort(const std::string& msg);
void print_configuration(std::ostream& os);
void finalize();
void finalize_kokkos();
void finalize_session();

}
