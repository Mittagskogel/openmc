#include "openmc/dagmc.h"

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/material.h"
#include "openmc/string_utils.h"
#include "openmc/settings.h"
#include "openmc/surface.h"

#ifdef DAGMC
#include "uwuw.hpp"
#include "dagmcmetadata.hpp"
#endif
#include <fmt/core.h>

#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>

namespace openmc {

#ifdef DAGMC
const bool DAGMC_ENABLED = true;
#else
const bool DAGMC_ENABLED = false;
#endif

}

#ifdef DAGMC

namespace openmc {

const std::string DAGMC_FILENAME = "dagmc.h5m";

std::string dagmc_file() {
  std::string filename = settings::path_input + DAGMC_FILENAME;
  if (!file_exists(filename)) {
    fatal_error("Geometry DAGMC file '" + filename + "' does not exist!");
  }
  return filename;
}

std::string
DAGUniverse::dagmc_ids_for_dim(int dim) const
{
  // generate a vector of ids
  std::vector<int> id_vec;
  int n_cells = dagmc_instance_->num_entities(dim);
  for (int i = 1; i <= n_cells; i++) {
    id_vec.push_back(dagmc_instance_->id_by_index(dim, i));
  }

  // sort the vector of ids
  std::sort(id_vec.begin(), id_vec.end());

  // generate a string representation of the range
  std::stringstream out;

  int i = 0;
  int start_id = id_vec[0];
  int stop_id;
  while (i < n_cells) {
    stop_id = id_vec[i];

    if (id_vec[i + 1] > stop_id + 1) {
      if (start_id != stop_id) {
        // there are several IDs in a row, print condensed version
        out << start_id << "-" << stop_id;
      } else {
        // only one ID in this contiguous block
        out << start_id;
      }
      if (i < n_cells - 1) { out << ", "; }
      start_id = id_vec[++i];
      stop_id = start_id;
    }

    i++;
  }

  return out.str();
}

bool
DAGUniverse::uses_uwuw() const
{
  return uwuw_ && !uwuw_->material_library.empty();
}

std::string
DAGUniverse::get_uwuw_materials_xml() const
{

  if (!uses_uwuw()) {
    throw std::runtime_error("This DAGMC Universe does not use UWUW materials");
  }

  std::stringstream ss;
  // write header
  ss << "<?xml version=\"1.0\"?>\n";
  ss << "<materials>\n";
  const auto& mat_lib = uwuw_->material_library;
  // write materials
  for (auto mat : mat_lib) { ss << mat.second->openmc("atom"); }
  // write footer
  ss << "</materials>";

  return ss.str();
}

void
DAGUniverse::write_uwuw_materials_xml(const std::string& outfile) const
{
  if (!uses_uwuw()) {
    throw std::runtime_error("This DAGMC universe does not use UWUW materials.");
  }

  std::string xml_str = get_uwuw_materials_xml();
  // if there is a material library in the file
  std::ofstream mats_xml(outfile);
  mats_xml << xml_str;
  mats_xml.close();
}

void
DAGUniverse::legacy_assign_material(std::string mat_string,
                                    std::unique_ptr<DAGCell>& c) const
{
  bool mat_found_by_name = false;
  // attempt to find a material with a matching name
  to_lower(mat_string);
  for (const auto& m : model::materials) {
    std::string m_name = m->name();
    to_lower(m_name);
    if (mat_string == m_name) {
      // assign the material with that name
      if (!mat_found_by_name) {
        mat_found_by_name = true;
        c->material_.push_back(m->id_);
      // report error if more than one material is found
      } else {
        fatal_error(fmt::format(
          "More than one material found with name '{}'. Please ensure materials "
          "have unique names if using this property to assign materials.",
          mat_string));
      }
    }
  }

  // if no material was set using a name, assign by id
  if (!mat_found_by_name) {
    try {
      auto id = std::stoi(mat_string);
      c->material_.emplace_back(id);
    } catch (const std::invalid_argument&) {
      fatal_error(fmt::format(
        "No material '{}' found for volume (cell) {}", mat_string, c->id_));
    }
  }

  if (settings::verbosity >= 10) {
    const auto& m = model::materials[model::material_map.at(c->material_[0])];
    std::stringstream msg;
    msg << "DAGMC material " << mat_string << " was assigned";
    if (mat_found_by_name) {
      msg << " using material name: " << m->name_;
    } else {
      msg << " using material id: " << m->id_;
    }
    write_message(msg.str(), 10);
  }
}

void read_dagmc_universes(pugi::xml_node node) {
  for (pugi::xml_node dag_node : node.children("dagmc")) {
    model::universes.push_back(std::make_unique<DAGUniverse>(dag_node));
    model::universe_map[model::universes.back()->id_] = model::universes.size() - 1;
  }
}

DAGUniverse::DAGUniverse(pugi::xml_node node) {
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify the id of the DAGMC universe");
  }

  if (check_for_node(node, "filename")) {
    filename_ = get_node_value(node, "filename");
  } else {
    fatal_error("Must specify a file for the DAGMC universe");
  }

  adjust_geometry_ids_ = false;
  if (check_for_node(node, "auto_geom_ids")) {
    adjust_geometry_ids_ = get_node_value_bool(node, "auto_geom_ids");
  }

  adjust_material_ids_ = false;
  if (check_for_node(node, "auto_mat_ids")) {
    adjust_geometry_ids_ = get_node_value_bool(node, "auto_mat_ids");
  }

  initialize();
}

DAGUniverse::DAGUniverse(const std::string& filename, bool auto_geom_ids)
: filename_(filename), adjust_geometry_ids_(auto_geom_ids) {
  // determine the next universe id
  int32_t next_univ_id = 0;
  for (const auto& u : model::universes) {
    if (u->id_ > next_univ_id) next_univ_id = u->id_;
  }
  next_univ_id++;

  // set the universe id
  id_ = next_univ_id;

  initialize();
}

void DAGUniverse::initialize() {
  type_ = GeometryType::DAG;

  // determine the next cell id
  int32_t next_cell_id = 0;
  for (const auto& c : model::cells) {
    if (c->id_ > next_cell_id) next_cell_id = c->id_;
  }
  cell_idx_offset_ = model::cells.size();
  next_cell_id++;

  // determine the next surface id
  int32_t next_surf_id = 0;
  for (const auto& s : model::surfaces) {
    if (s->id_ > next_surf_id) next_surf_id = s->id_;
  }
  surf_idx_offset_ = model::surfaces.size();
  next_surf_id++;

  // create a new DAGMC instance
  dagmc_instance_ = std::make_shared<moab::DagMC>();

  // --- Materials ---

  // read any UWUW materials from the file
  read_uwuw_materials();

  // check for uwuw material definitions
  bool using_uwuw = uses_uwuw();

  // notify user if UWUW materials are going to be used
  if (using_uwuw) {
    write_message("Found UWUW Materials in the DAGMC geometry file.", 6);
  }

  // load the DAGMC geometry
  moab::ErrorCode rval = dagmc_instance_->load_file(filename_.c_str());
  MB_CHK_ERR_CONT(rval);

  // initialize acceleration data structures
  rval = dagmc_instance_->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  // parse model metadata
  dagmcMetaData DMD(dagmc_instance_.get(), false, false);
  DMD.load_property_data();

  std::vector<std::string> keywords {"temp"};
  std::map<std::string, std::string> dum;
  std::string delimiters = ":/";
  rval = dagmc_instance_->parse_properties(keywords, dum, delimiters.c_str());
  MB_CHK_ERR_CONT(rval);

  // --- Cells (Volumes) ---

  // initialize cell objects
  int n_cells = dagmc_instance_->num_entities(3);
  moab::EntityHandle graveyard = 0;
  for (int i = 0; i < n_cells; i++) {
    moab::EntityHandle vol_handle = dagmc_instance_->entity_by_index(3, i + 1);

    // set cell ids using global IDs
    auto c = std::make_unique<DAGCell>();
    c->dag_index_ = i + 1;
    c->id_ = adjust_geometry_ids_ ? next_cell_id++ : dagmc_instance_->id_by_index(3, c->dag_index_);
    c->dagmc_ptr_ = dagmc_instance_;
    c->universe_ = id_;
    c->fill_ = C_NONE; // no fill, single universe

   auto in_map = model::cell_map.find(c->id_);
    if (in_map == model::cell_map.end()) {
      model::cell_map[c->id_] = model::cells.size();
    } else {
      warning(fmt::format("DAGMC Cell IDs: {}", dagmc_ids_for_dim(3)));
      fatal_error(fmt::format("Cell ID {} exists in both DAGMC Universe {} "
                              "and the CSG geometry.", c->id_, this->id_));
    }

    // MATERIALS

    // determine volume material assignment
    std::string mat_str = DMD.get_volume_property("material", vol_handle);

    if (mat_str.empty()) {
      fatal_error(fmt::format("Volume {} has no material assignment.", c->id_));
    }

    to_lower(mat_str);

    if (mat_str == "graveyard") {
      graveyard = vol_handle;
    }

    // material void checks
    if (mat_str == "void" || mat_str == "vacuum" || mat_str == "graveyard") {
      c->material_.push_back(MATERIAL_VOID);
    } else {
      if (using_uwuw) {
        // lookup material in uwuw if present
        std::string uwuw_mat = DMD.volume_material_property_data_eh[vol_handle];
        if (uwuw_->material_library.count(uwuw_mat) != 0) {
          // Note: material numbers are set by UWUW
          int mat_number = uwuw_->material_library.get_material(uwuw_mat).metadata["mat_number"].asInt();
          c->material_.push_back(mat_number);
        } else {
          fatal_error(fmt::format("Material with value '{}' not found in the "
            "UWUW material library", mat_str));
        }
      } else {
        legacy_assign_material(mat_str, c);
      }
    }

    // check for temperature assignment
    std::string temp_value;

    // no temperature if void
    if (c->material_[0] == MATERIAL_VOID) {
      model::cells.emplace_back(std::move(c));
      continue;
    }

    // assign cell temperature
    const auto& mat = model::materials[model::material_map.at(c->material_[0])];
    if (dagmc_instance_->has_prop(vol_handle, "temp")) {
      rval = dagmc_instance_->prop_value(vol_handle, "temp", temp_value);
      MB_CHK_ERR_CONT(rval);
      double temp = std::stod(temp_value);
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * temp));
    } else if (mat->temperature() > 0.0) {
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * mat->temperature()));
    } else {
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * settings::temperature_default));
    }

    model::cells.emplace_back(std::move(c));
  }

  // allocate the cell overlap count if necessary
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

  if (!graveyard) {
    warning("No graveyard volume found in the DagMC model. "
            "This may result in lost particles and rapid simulation failure.");
  }

  // --- Surfaces ---

  // initialize surface objects
  int n_surfaces = dagmc_instance_->num_entities(2);
  for (int i = 0; i < n_surfaces; i++) {
    moab::EntityHandle surf_handle = dagmc_instance_->entity_by_index(2, i+1);

    // set cell ids using global IDs
    auto s = std::make_unique<DAGSurface>();
    s->dag_index_ = i+1;
    s->id_ = adjust_geometry_ids_ ? next_surf_id++ : dagmc_instance_->id_by_index(2, i+1);
    s->dagmc_ptr_ = dagmc_instance_;

    // set BCs
    std::string bc_value = DMD.get_surface_property("boundary", surf_handle);
    to_lower(bc_value);
    if (bc_value.empty() || bc_value == "transmit" || bc_value == "transmission") {
      // set to transmission by default (nullptr)
    } else if (bc_value == "vacuum") {
      s->bc_ = std::make_shared<VacuumBC>();
    } else if (bc_value == "reflective" || bc_value == "reflect" || bc_value == "reflecting") {
      s->bc_ = std::make_shared<ReflectiveBC>();
    } else if (bc_value == "periodic") {
      fatal_error("Periodic boundary condition not supported in DAGMC.");
    } else {
      fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
        "on surface {}", bc_value, s->id_));
    }

    // graveyard check
    moab::Range parent_vols;
    rval = dagmc_instance_->moab_instance()->get_parent_meshsets(surf_handle, parent_vols);
    MB_CHK_ERR_CONT(rval);

    // if this surface belongs to the graveyard
    if (graveyard && parent_vols.find(graveyard) != parent_vols.end()) {
      // set graveyard surface BC's to vacuum
      s->bc_ = std::make_shared<VacuumBC>();
    }

    // add to global array and map

    auto in_map = model::surface_map.find(s->id_);
    if (in_map == model::surface_map.end()) {
      model::surface_map[s->id_] = model::surfaces.size();
    } else {
      warning(fmt::format("DAGMC Surface IDs: {}", dagmc_ids_for_dim(2)));
      fatal_error(fmt::format("Surface ID {} exists in both Universe {} "
                              "and the CSG geometry.", s->id_, this->id_));
    }

    model::surfaces.emplace_back(std::move(s));

  } // end surface loop
}

void
DAGUniverse::read_uwuw_materials() {

  int32_t next_material_id = 0;
  for (const auto& m : model::materials) {
    next_material_id = std::max(m->id_, next_material_id);
  }
  next_material_id++;

  uwuw_ = std::make_shared<UWUW>(filename_.c_str());
  const auto& mat_lib = uwuw_->material_library;
  if (mat_lib.size() == 0) return;

  // if we're using automatic IDs, update the UWUW material metadata
  if (adjust_material_ids_) {
    for (auto& mat : uwuw_->material_library) {
      mat.second->metadata["mat_number"] = std::to_string(next_material_id++);
    }
  }

  std::stringstream ss;
  ss << "<?xml version=\"1.0\"?>\n";
  ss << "<materials>\n";
  for (auto mat : mat_lib) { ss << mat.second->openmc("atom"); }
  ss << "</materials>";
  std::string mat_xml_string = ss.str();

  // create a pugi XML document from this string
  pugi::xml_document doc;
  auto result = doc.load_string(mat_xml_string.c_str());
  if (!result) {
    fatal_error("Error processing XML created using DAGMC UWUW materials.");
  }
  pugi::xml_node root = doc.document_element();
  for (pugi::xml_node material_node : root.children("material")) {
    model::materials.push_back(std::make_unique<Material>(material_node));
  }
}

}
#endif
