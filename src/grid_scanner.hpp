// ====================================================================
// Provides a wrapper for doing grid scans with an arbitrary
// number of parameters. Mainly just for convenience, since this
// gets done so often.
// ====================================================================

#ifndef GRID_SCANNER_H
#define GRID_SCANNER_H

#include <cassert>
#include <vector>

namespace flexiblesusy {

/**
 * @class Grid_scanner
 * @brief Does grid scan over an arbitrary dimension space
 *
 * Counting is zero-based, i.e. the first point has each
 * index 0. Periodic boundary conditions are imposed on
 * the linearised grid, so that stepping forward past
 * the end of the linear array wraps around to the
 * first grid point. Likewise for the grid itself,
 * so that requesting a position at an index greater
 * than the number of points in that direction wraps
 * around. Note that once has_finished() returns true
 * the current position will be at the start of the
 * grid once again (i.e. no need to step forward an
 * extra step).
 *
 * Example:
 * @code
 * struct MyExample {
 *    static const std::vector<double> func(std::vector<std::size_t> grid_pt)
 *       {
 *         double x_incr = 1.5;
 *         double y_incr = 2.0;
 *         double z_incr = 3.0;
 *         std::vector<double> temp(grid_pt.size(), 0.0);
 *        
 *         temp.at(0) = 2.0 + x_incr * grid_pt.at(0);
 *         temp.at(1) = 1.0 + y_incr * grid_pt.at(1);
 *         temp.at(2) = z_incr * grid_pt.at(2);
 * 
 *         return temp; 
 *       }
 * };
 * 
 *
 * std::vector<std::size_t> dims = {3, 2, 5};
 * Grid_scanner scan(dims);
 * std::vector<double> posn;
 * while (!scan.has_finished()) {
 *    posn = MyExample::func(scan.get_position());
 *    std::cout << "x(" << scan.get_current_index() << ") = (";
 *    for (std::size_t i = 0; i < posn.size(); ++i) {
 *       std::cout << posn.at(i);
 *       if (i != posn.size() - 1) {
 *          std::cout << ", ";
 *       } else {
 *          std::cout << ")\n";
 *       }
 *    }
 *    scan.step_forward();
 * }
 * 
 * @endcode
 */
   class Grid_scanner {
   public:
      
      Grid_scanner(const std::vector<std::size_t>&);

      std::size_t get_num_dims() const { return num_dims; }
      std::size_t get_num_points() const { return num_points; }
      std::size_t get_current_index() const;
      bool has_finished() const { return finished; };
      const std::vector<std::size_t>& get_position() const
         {
            return position;
         }    
      const std::vector<std::size_t>& get_dimensions() const
         {
            return dimensions;
         }

      void set_position(const std::vector<std::size_t>&);

      void reset();
      void step_forward();
      void step_backward();

   private:
      std::size_t num_dims; ///< number of dimensions in the grid
      std::size_t num_points; ///< number of grid points
      bool finished; ///< flag indicating have reached end of scan
      std::vector<std::size_t> position; ///< current position in the grid
      std::vector<std::size_t> dimensions; ///< number of points along each direction
      std::vector<std::size_t> coefficients; ///< coefficients in conversion to linear index

   };

} // namespace flexiblesusy

#endif
