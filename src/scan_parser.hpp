// ====================================================================
// Factored out parsing code for doing scans over input parameters
// ====================================================================

#ifndef SCAN_PARSER_H
#define SCAN_PARSER_H

#include <map>
#include <string>
#include <vector>

namespace flexiblesusy {

struct Scanned_parameter {
   std::string name;
   double lower_bound;
   double upper_bound;
   std::size_t num_points;

   Scanned_parameter(const std::string& name_)
      : name(name_), lower_bound(0.), upper_bound(0.)
      , num_points(1)
      {
      }
};

struct Fixed_parameter {
   std::string name;
   double value;

   Fixed_parameter(const std::string& name_)
      : name(name_), value(0.)
      {
      }
};

class Scan_parser {
public:
   Scan_parser();
   ~Scan_parser();

   std::vector<Scanned_parameter> get_scanned_parameters() const;
   std::vector<Fixed_parameter> get_fixed_parameters() const;
   bool is_grid_scan() const { return grid_scan; }
   int get_number_of_points() const { return total_npts; }

   void parse_scan_inputs_file(const std::string&);

private:
   bool grid_scan;
   std::map<std::string,Scanned_parameter> scanned_parameters;
   std::map<std::string,Fixed_parameter> fixed_parameters;
   int total_npts;

   std::string lower_bnd_suffix;
   std::string upper_bnd_suffix;
   std::string num_points_suffix;
   std::string value_suffix;

   std::size_t calculate_total_num_points() const;
   void check_scan_bounds();
   bool ends_with(const std::string&, const std::string&) const;
   void set_fixed_parameter_property(const std::string&, const std::string&, const std::string&);
   void set_scanned_parameter_property(const std::string&, const std::string&, const std::string&);
   void trim(std::string&) const;
};

} // namespace flexiblesusy

#endif
