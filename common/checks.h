#pragma once

#include <iostream>
#include <stdexcept>

// Check for nullptr
#define __check_not_null_ret(var, msg) if (var == nullptr) { std::cerr << msg << std::endl; return 0; }
#define __check_not_null_throw(var, msg) if (var == nullptr) { throw std::runtime_error(msg); }

// Check number of components
#define __check_num_components_ret(var, num_comp, msg) if (var->GetNumberOfComponents() != num_comp) { std::cerr << msg << std::endl; return 0; }
#define __check_num_components_throw(var, num_comp, msg) if (var->GetNumberOfComponents() != num_comp) { throw std::runtime_error(msg); }

// Check number of stored values
#define __check_not_empty_ret(var, msg, ok) if (var->GetNumberOfValues() == 0) { if (ok) { std::cout << msg << std::endl; } else { std::cerr << msg << std::endl; } return ok; }
#define __check_not_empty_throw(var, msg) if (var->GetNumberOfValues() == 0) { throw std::runtime_error(msg); }

// Check parameters for minimum and maximum values
#define __check_param_min_ret(var, min_val, msg, ok) if (var < min_val) { if (ok) { std::cout << msg << std::endl; } else { std::cerr << msg << std::endl; } return ok; }
#define __check_param_min_throw(var, min_val, msg) if (var < min_val) { throw std::runtime_error(msg); }

#define __check_param_max_ret(var, max_val, msg, ok) if (var > max_val) { if (ok) { std::cout << msg << std::endl; } else { std::cerr << msg << std::endl; } return ok; }
#define __check_param_max_throw(var, max_val, msg) if (var > max_val) { throw std::runtime_error(msg); }
