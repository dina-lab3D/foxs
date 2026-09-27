#include "CIF.h"

#include <boost/algorithm/string.hpp>
#include <iostream>

std::map<std::string, std::size_t> CIF::atomSiteTable_;

bool CIF::readAtomSiteTable(std::istream& ifile)
{
  unsigned int index = 0;
  int place = 0;
  atomSiteTable_.clear();
  while (!ifile.eof()) {
    std::string line;
    getline(ifile, line);
    if(isAtomSiteLine(line)) {
      boost::trim(line);
      std::string key = line.substr(11);
      atomSiteTable_[key] = index;
      index++;
    } else {
      // we are done reading the atom site dictionary, go back to the beg. of the line
      if(index > 0) {
        ifile.seekg(place);
        break;
      }
    }
    place = ifile.tellg();
  }
  if(index > 0) return true;
  return false;
}

void CIF::writeAtomSiteTable(std::ostream& ofile) {
  ofile << "#\n loop_\n";

  ofile << "_atom_site.group_PDB " << std::endl;
  ofile << "_atom_site.id " << std::endl;
  ofile << "_atom_site.type_symbol " << std::endl;
  ofile << "_atom_site.label_atom_id " << std::endl;
  ofile << "_atom_site.label_alt_id " << std::endl;
  ofile << "_atom_site.label_comp_id " << std::endl;
  ofile << "_atom_site.label_asym_id " << std::endl;
  ofile << "_atom_site.label_entity_id " << std::endl;
  ofile << "_atom_site.label_seq_id " << std::endl;
  ofile << "_atom_site.pdbx_PDB_ins_code " << std::endl;
  ofile << "_atom_site.Cartn_x " << std::endl;
  ofile << "_atom_site.Cartn_y " << std::endl;
  ofile << "_atom_site.Cartn_z " << std::endl;
  ofile << "_atom_site.occupancy " << std::endl;
  ofile << "_atom_site.B_iso_or_equiv " << std::endl;
  ofile << "_atom_site.pdbx_formal_charge " << std::endl;
  ofile << "_atom_site.auth_seq_id " << std::endl;
  ofile << "_atom_site.auth_comp_id " << std::endl;
  ofile << "_atom_site.auth_asym_id " << std::endl;
  ofile << "_atom_site.auth_atom_id " << std::endl;
  ofile << "_atom_site.pdbx_PDB_model_num " << std::endl;
}

bool CIF::isAtomSiteLine(const std::string& cifLine)
{
  return cifLine.substr(0,11) == "_atom_site.";
}

std::string CIF::getColumn(const std::string& cifLine, const std::string& key)
{
  std::vector<std::string> result;
  boost::split(result, cifLine, boost::is_any_of("\t "), boost::token_compress_on);
  if(atomSiteTable_.find(key) != atomSiteTable_.end()) {
    unsigned int col = atomSiteTable_[key];
    if(col < result.size())
      return result[col];
  }
  return std::string();
}

bool CIF::isATOMrec(const std::string& cifLine)
{
  return getColumn(cifLine, "group_PDB") == "ATOM";
}

bool CIF::isHETATMrec(const std::string& cifLine)
{
  return getColumn(cifLine, "group_PDB") == "HETATM";
}

float CIF::atomXCoord(const std::string& cifLine)
{
  try {
    return std::stof(getColumn(cifLine, "Cartn_x"));
  } catch(...) {
    return 0.0;
  }
}

float CIF::atomYCoord(const std::string& cifLine)
{
  try {
    return std::stof(getColumn(cifLine, "Cartn_y"));
  } catch(...) {
    return 0.0;
  }
}

float CIF::atomZCoord(const std::string& cifLine)
{
  try {
    return std::stof(getColumn(cifLine, "Cartn_z"));
  } catch(...) {
    return 0.0;
  }
}

int CIF::atomIndex(const std::string& cifLine)
{
  try {
    return std::stof(getColumn(cifLine, "id"));
  } catch(...) {
    return 0;
  }
}

std::string CIF::atomChainId(const std::string& cifLine)
{
  return getColumn(cifLine, "label_asym_id");
}

std::string CIF::atomChainIdAuthor(const std::string& cifLine)
{
  return getColumn(cifLine, "auth_asym_id");
}

int CIF::atomResidueIndex(const std::string& cifLine)
{
  int ret = 0;
  try {
    ret = std::stoi(getColumn(cifLine, "auth_seq_id"));
  } catch(std::exception& err) {
    try {
      ret = std::stoi(getColumn(cifLine, "label_seq_id"));
    }
    catch(std::exception& err) {
      return ret;
    }
  }
  return ret;
}

std::string CIF::atomType(const std::string& cifLine)
{
  return getColumn(cifLine, "label_atom_id");
}

std::string CIF::getAtomResidueLongName(const std::string& cifLine)
{
  return getColumn(cifLine, "label_comp_id");
}

AtomEntryType CIF::getAtomEntryType(const std::string& cifLine)
{
  if(isHETATMrec(cifLine))
    return AtomEntryType::HETATM;
  else if (isATOMrec(cifLine))
    return AtomEntryType::ATOM;
  else
    return AtomEntryType::UNK;
}

float CIF::getOccupancy(const std::string& cifLine)
{
  return std::stof(getColumn(cifLine,"occupancy"));
}

float CIF::getTempFactor(const std::string& cifLine)
{
  return std::stof(getColumn(cifLine,"B_iso_or_equiv"));
}

int CIF::getModelNum(const std::string& cifLine) {
  return std::stoi(getColumn(cifLine,"pdbx_PDB_model_num"));
}
