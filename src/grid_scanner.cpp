#include "grid_scanner.hpp"

namespace flexiblesusy {

   Grid_scanner::Grid_scanner(const std::vector<std::size_t>& dims)
      : num_dims(dims.size())
      , finished(false)
      , position(dims.size(), 0)
      , dimensions(dims)
   {
      assert(dims.size() >= 1 && "Grid_scanner::Grid_scanner: vector must have at"
             " least one dimension!");

      num_points = 1;
      coefficients = std::vector<std::size_t>(dims.size(), 1);
      for (std::size_t i = dims.size() - 1; i > 0; --i) {
         assert(dims[i] > 0 && "Grid_scanner::Grid_scanner: number of points in dimension"
                " must be at least one!");

         num_points *= dims[i];
         coefficients[i-1] = dims[i] * coefficients[i];
      }

      assert(dims[0] > 0 && "Grid_scanner::Grid_scanner: number of points in dimension"
             " must be at least one!");
      
      num_points *= dims[0];
   }

   std::size_t Grid_scanner::get_current_index() const
   {
      std::size_t index = 0;

      for (std::size_t i = 0; i < num_dims; ++i) {
         index += position[i] * coefficients[i];
      }

      return index;
   }

   void Grid_scanner::set_position(const std::vector<std::size_t>& posn)
   {
      assert(posn.size() == num_dims && "Grid_scanner::set_position: length of vector must match"
             " number of dimensions!");

      for (std::size_t i = 0; i < num_dims; ++i) {
         position[i] = posn[i] % dimensions[i];
      }
   }

   void Grid_scanner::reset() 
   {
      for (std::size_t i = 0; i < num_dims; ++i) {
         position[i] = 0;
      }
      finished = false;
   }

   void Grid_scanner::step_forward()
   {
      finished = true;

      for (std::size_t i = num_dims; i > 0; --i) {
         ++position[i-1];
         if (position[i-1] != dimensions[i-1]) {
            finished = false;
            break;
         }
         position[i-1] = 0;
      }
      
   }

   void Grid_scanner::step_backward()
   {
      finished = true;

      for (std::size_t i = num_dims; i > 0; --i) {
         if (position[i-1] != 0) {
            --position[i-1];
            finished = false;
            break;
         }
         position[i-1] = dimensions[i-1] - 1;
      }
   }

} // namespace flexiblesusy
