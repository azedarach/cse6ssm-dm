// ====================================================================
// Implementation of scan parser, factored out from old scan code
// ====================================================================

#include "scan_parser.hpp"

#include "error.hpp"
#include "logger.hpp"

#include <fstream>
#include <limits>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace flexiblesusy {

Scan_parser::Scan_parser()
   : grid_scan(true)
   , scanned_parameters()
   , fixed_parameters()
   , total_npts(1)
   , lower_bnd_suffix("_lower")
   , upper_bnd_suffix("_upper")
   , num_points_suffix("_npts")
   , value_suffix("_value")
{
}

Scan_parser::~Scan_parser()
{
}

std::size_t Scan_parser::calculate_total_num_points() const
{
   std::size_t num_points = 1;
   for (std::map<std::string,Scanned_parameter>::const_iterator it = scanned_parameters.begin(),
           end = scanned_parameters.end(); it != end; ++it) {
      num_points *= it->second.num_points;
   }

   return num_points;
}

void Scan_parser::check_scan_bounds()
{
   for (std::map<std::string,Scanned_parameter>::iterator it = scanned_parameters.begin(),
           end = scanned_parameters.end(); it != end; ++it) {
      if (it->second.lower_bound > it->second.upper_bound) {
         double tmp = it->second.lower_bound;
         it->second.lower_bound = it->second.upper_bound;
         it->second.upper_bound = tmp;
      }
   }
}

bool Scan_parser::ends_with(const std::string& str, const std::string& suffix) const
{
   return !str.compare(str.size() - suffix.size(), std::string::npos, suffix);
}

std::vector<Fixed_parameter> Scan_parser::get_fixed_parameters() const
{
   std::vector<Fixed_parameter> params;
   for (std::map<std::string,Fixed_parameter>::const_iterator it = fixed_parameters.begin(),
           end = fixed_parameters.end(); it != end; ++it) {
      params.push_back(it->second);
   }
   return params;
}

std::vector<Scanned_parameter> Scan_parser::get_scanned_parameters() const
{
   std::vector<Scanned_parameter> params;
   for (std::map<std::string,Scanned_parameter>::const_iterator it = scanned_parameters.begin(),
           end = scanned_parameters.end(); it != end; ++it) {
      params.push_back(it->second);
   }
   return params;
}

void Scan_parser::parse_scan_inputs_file(const std::string& input_file)
{
   std::ifstream ifs(input_file, std::ifstream::in);
   if (ifs.fail())
      throw ReadError("unable to open file " + input_file);

   bool total_npts_requested = false;

   std::string line;
   while (std::getline(ifs, line)) {
      if (line.find("#") != std::string::npos)
         line = line.substr(0, line.find("#"));
      if (line.empty())
         continue;

      // break up into individual fields
      std::vector<std::string> fields;
      boost::split(fields, line, boost::is_any_of(";"));

      for (std::size_t i = 0; i < fields.size(); ++i) {
         // remove whitespace
         trim(fields[i]);

         // get field description and value
         if (!fields[i].empty()) {
            std::vector<std::string> field;
            boost::split(field, fields[i], boost::is_any_of("="));

            if (field.size() != 2) {
               WARNING("Ignoring invalid input '" + fields[i] + "'");
               continue;
            }
             
            trim(field[0]);
            trim(field[1]);
            
            // compare against valid inputs
            if (field[0] == "is_grid_scan") {
               boost::to_lower(field[1]);
               if (field[1] == "false" || field[1] == "f") {
                  grid_scan = false;
               } else if (field[1] == "true" || field[1] == "t") {
                  grid_scan = true;
               } else {
                  WARNING("Ignoring invalid input '" + fields[i] + "'");
               }
            } else if (field[0] == "total_npts") {
               int num_points;
               try {
                  num_points = boost::lexical_cast<int>(field[1]);
                  total_npts_requested = true;
                  if (num_points <= 0) {
                     WARNING("Ignoring invalid input '" + fields[i] + "'");
                     num_points = 1;
                     total_npts_requested = false;
                  }
                  total_npts = num_points;
               } catch (const boost::bad_lexical_cast& error) {
                  WARNING("Ignoring invalid input '" + fields[i] + "'");
                  total_npts = 1;
                  total_npts_requested = false;
               }
            } else if (ends_with(field[0], lower_bnd_suffix)) {
               set_scanned_parameter_property(field[0], field[1], lower_bnd_suffix);
            } else if (ends_with(field[0], upper_bnd_suffix)) {
               set_scanned_parameter_property(field[0], field[1], upper_bnd_suffix);               
            } else if (ends_with(field[0], num_points_suffix)) {
               set_scanned_parameter_property(field[0], field[1], num_points_suffix);
            } else if (ends_with(field[0], value_suffix)) {
               set_fixed_parameter_property(field[0], field[1], value_suffix);
            } else {
               WARNING("Ignoring invalid input '" + fields[i] + "'");
            }
         }
      }
   }
   
   if (total_npts_requested && grid_scan) {
      WARNING("Grid scan requested, total number of points will be ignored\n");
   }

   if (grid_scan) {
      total_npts = calculate_total_num_points();
   }

   check_scan_bounds();
}

void Scan_parser::set_fixed_parameter_property(const std::string& desc, const std::string& value, const std::string& prop)
{
   const std::string name = desc.substr(0, desc.size() - prop.size());

   Fixed_parameter tmp_param(name);
   tmp_param.value = 0.;

   std::pair<std::map<std::string,Fixed_parameter>::iterator,bool> param
      = fixed_parameters.insert(std::pair<std::string,Fixed_parameter>(name, tmp_param));

   if (prop == value_suffix) {
      double val = 0.;
      try {
         val = boost::lexical_cast<double>(value);
      } catch (const boost::bad_lexical_cast& error) {
         WARNING("Ignoring invalid lower bound '" + value + "' for field '" + desc + "'");
         val = 0.;
      }
      param.first->second.value = val;
   } else {
      WARNING("Unrecognised property requested!");
   }
}

void Scan_parser::set_scanned_parameter_property(const std::string& desc, const std::string& value, const std::string& prop)
{
   const std::string name = desc.substr(0, desc.size() - prop.size());

   Scanned_parameter tmp_param(name);

   // by default, if no lower or upper bounds are given
   // the parameter is unbounded
   tmp_param.lower_bound = -std::numeric_limits<double>::max();
   tmp_param.upper_bound = std::numeric_limits<double>::max();
   tmp_param.num_points = 1;

   std::pair<std::map<std::string,Scanned_parameter>::iterator,bool> param
      = scanned_parameters.insert(std::pair<std::string,Scanned_parameter>(name, tmp_param));

   if (prop == lower_bnd_suffix) {
      double lower_bound = -std::numeric_limits<double>::max();
      try {
         lower_bound = boost::lexical_cast<double>(value);
      } catch (const boost::bad_lexical_cast& error) {
         WARNING("Ignoring invalid lower bound '" + value + "' for field '" + desc + "'");
         lower_bound = -std::numeric_limits<double>::max();
      }
      param.first->second.lower_bound = lower_bound;
   } else if (prop == upper_bnd_suffix) {
      double upper_bound = std::numeric_limits<double>::max();
      try {
         upper_bound = boost::lexical_cast<double>(value);
      } catch (const boost::bad_lexical_cast& error) {
         WARNING("Ignoring invalid upper bound '" + value + "' for field '" + desc + "'");
         upper_bound = std::numeric_limits<double>::max();
      }
      param.first->second.upper_bound = upper_bound;
   } else if (prop == num_points_suffix) {
      int num_points;
      try {
         num_points = boost::lexical_cast<int>(value);
         if (num_points <= 0) {
            WARNING("Ignoring invalid number of points '" + value + "' for field '" + desc + "'");
            num_points = 1;
         }
      } catch (const boost::bad_lexical_cast& error) {
         WARNING("Ignoring invalid number of points '" + value + "' for field '" + desc + "'");
         num_points = 1;
      }
      param.first->second.num_points = num_points;
   } else {
      WARNING("Unrecognised property requested!");
   }
}

void Scan_parser::trim(std::string& str) const
{
   std::size_t startpos = str.find_first_not_of(" \t\n\v\f\r");
   if (startpos != std::string::npos) str.erase(0, startpos);
   std::size_t endpos = str.find_last_not_of(" \t\n\v\f\r");
   if (endpos != std::string::npos) str.erase(endpos+1);
}

} // namespace flexiblesusy
