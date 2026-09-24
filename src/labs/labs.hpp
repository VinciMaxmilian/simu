#pragma once
#include <memory>
#include <string>

#include "core/lab.hpp"

std::unique_ptr<Lab> make_toy_lab(uint64_t seed);
std::unique_ptr<Lab> make_chemistry_lab();
std::unique_ptr<Lab> make_materials_lab();
std::unique_ptr<Lab> make_nuclear_lab();
