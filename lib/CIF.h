#ifndef _CIF_h
#define _CIF_h

#include <string>
#include <iosfwd>
#include <map>

#include "PDB.h"

/*
CLASS
  CIF

  The CIF class defines a set of tools used for reading the Protein Data Bank
  file format, used to describe molecular structures, like proteins,
  drugs and Nucleic Acids

  KEYWORDS
  CIF, protein, Nucleic Acids (RNA, DNA), atom, residue, structure, coordinate, field

AUTHORS
  Dina Schneidman

GOALS
  Class CIF is designed to encapsulate a set of methods used for reading the
  CIF file format. Methods are static requiring no object instantiation
  before their use.

USAGE
  Similar to PDB
  Selectors: uses PDB selectors by converting CIF line to PDB line
*/
class CIF {

public:

  // fill in the atom_site columns table/dictionary
  static bool readAtomSiteTable(std::istream& ifile);
  static void writeAtomSiteTable(std::ostream& ofile);

  //
  static bool isAtomSiteLine(const std::string& cifLine);

  // Returns true if the given record is an ATOM record.
  static bool isATOMrec(const std::string& cifLine);

  // Returns true if the given record is a HETATM record.
  static bool isHETATMrec(const std::string& cifLine);

  //// Returns an ATOM record's X coordinate.
  static float atomXCoord(const std::string& cifLine);

  //// Returns an ATOM record's Y coordinate.
  static float atomYCoord(const std::string& cifLine);

  //// Returns an ATOM record's Z coordinate.
  static float atomZCoord(const std::string& cifLine);

  //// Returns an ATOM record's atom index number.
  static int atomIndex(const std::string& cifLine);

  //// Returns atomType as a string
  static std::string atomType(const std::string& cifLine);

  //// Returns an ATOM record's chain id.
  static std::string atomChainId(const std::string& cifLine);

  //// Chain id named by author
  static std::string atomChainIdAuthor(const std::string& cifLine);

  //// Returns an ATOM record's residue index.
  static int atomResidueIndex(const std::string& cifLine);

  //// Returns the full string of the Residue Name field
  static std::string getAtomResidueLongName(const std::string& cifLine);

  //// Return the type of the CIF entry, currently only ATOM and HETATM are identified.
  static AtomEntryType getAtomEntryType(const std::string& cifLine);

  //// Returns the occupancy,  which is a measure of the fraction of molecules in the
  //   crystal in which the current atom actually occupies the specified position
  static float getOccupancy(const std::string& cifLine);

  //// Returns the temperature factor,  which is a measure of how much an atom
  // oscillates or vibrates around the specified position
  static float getTempFactor(const std::string& cifLine);

  static int getModelNum(const std::string& cifLine);

  static std::string getColumn(const std::string& cifLine, const std::string& key);

private:
  // Private constructor to prevent object instantiation. Class was meant to
  // be use as a static method pool.
  CIF() = default;

  // A CIF dictionary for ATOM records: key = field_name, value = column number
  static std::map<std::string, std::size_t> atomSiteTable_;

};

#endif
