#pragma once

// Safely delete object
#define __delete_ptr(ptr) if (ptr != nullptr) { delete ptr; ptr = nullptr; }
#define __delete_ptr_array(ptr) if (ptr != nullptr) { delete[] ptr; ptr = nullptr; }
#define __delete_ptr_vtk(ptr) if (ptr != nullptr) { ptr->Delete(); ptr = nullptr; }
