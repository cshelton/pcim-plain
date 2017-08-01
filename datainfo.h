#ifndef DATAINFO_H
#define DATAINFO_H

#include <string>
#include "serial.h"

// meta data about the data set (basically names of the labels/variables)
// useful for printing
struct datainfo {
	std::vector<std::string> svarnames, dvarnames;
	std::map<std::string,int> svarid,dvarid;
	inline std::string dvarname(int v) const {return v<0 ? "X" : dvarnames[v];}
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void save(Ar &ar, const unsigned int ver) const {
		ar & BOOST_SERIALIZATION_NVP(svarnames);
		ar & BOOST_SERIALIZATION_NVP(dvarnames);
	}
	template<typename Ar>
	void load(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(svarnames);
		ar & BOOST_SERIALIZATION_NVP(dvarnames);
		svarid.clear();
		for(int i=0;i<svarnames.size();i++) svarid[svarnames[i]] = i;
		dvarid.clear();
		for(int i=0;i<dvarnames.size();i++) dvarid[dvarnames[i]] = i;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};
BOOST_CLASS_EXPORT_KEY(datainfo)

#endif
