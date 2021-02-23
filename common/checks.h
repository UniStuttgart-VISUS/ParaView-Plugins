#pragma once

#include <iostream>
#include <stdexcept>

#define __check_not_null_ret(var, msg) if (var == nullptr) { std::cerr << msg << std::endl; return 0; }
#define __check_not_null_throw(var, msg) if (var == nullptr) { throw std::runtime_error(msg); }

#define __check_num_components_ret(var, num_comp, msg) if (var->GetNumberOfComponents() != num_comp) { std::cerr << msg << std::endl; return 0; }
#define __check_num_components_throw(var, num_comp, msg) if (var->GetNumberOfComponents() != num_comp) { throw std::runtime_error(msg); }
