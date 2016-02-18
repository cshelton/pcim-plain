#ifndef VARTRAJ_H
#define VARTRAJ_H

#include <map>
#include <boost/range/adaptor/reversed.hpp>
#include <vector>
#include <iostream>
#include <queue>
#include "serial.h"

class vartraj {
public:
	struct rec {
	public:
		rec() = default;
		rec(double val) { v = val; }
		double v;
		operator double() const { return v; }
	private:
		friend class boost::serialization::access;
		template<typename Ar>
		void serialize(Ar &ar, const unsigned int ver) {
			ar & BOOST_SERIALIZATION_NVP(v);
		}
	};
	typedef std::map<double, rec> trajT;
	typedef trajT::const_iterator const_iterator;

	//vartraj() = default;
	vartraj() : tstart(0.0), tend(0.0) {}

	vartraj(std::initializer_list<std::pair<double,double>> data)
	{for(auto & x : data) vtraj.insert(std::make_pair(x.first, rec(x.second)));}

	trajT::iterator begin() { return vtraj.begin(); }
	trajT::const_iterator begin() const { return vtraj.cbegin(); }
	trajT::const_iterator cbegin() const { return vtraj.cbegin(); }
	trajT::iterator end() { return vtraj.end(); }
	trajT::const_iterator end() const { return vtraj.cend(); }
	trajT::const_iterator cend() const { return vtraj.cend(); }

	trajT::reverse_iterator rbegin() { return vtraj.rbegin(); }
	trajT::const_reverse_iterator rbegin() const { return vtraj.crbegin(); }
	trajT::const_reverse_iterator crbegin() const { return vtraj.crbegin(); }
	trajT::reverse_iterator rend() { return vtraj.rend(); }
	trajT::const_reverse_iterator rend() const { return vtraj.crend(); }
	trajT::const_reverse_iterator crend() const { return vtraj.crend(); }

	std::size_t size() const {return vtraj.size();}
	bool empty() const {return vtraj.empty();}

	trajT::const_iterator find(double t) const {return vtraj.find(t);}
	trajT::const_iterator lower_bound(double t) const
	{return vtraj.lower_bound(t);}
	trajT::const_iterator upper_bound(double t) const
	{return vtraj.upper_bound(t);}

	vartraj(const vartraj & t) : vtraj(t.vtraj) {
		tstart = t.tstart;
		tend = t.tend;
	}

	vartraj(const vartraj & t, double endt) {
		tstart = t.tstart;
		tend = endt;
		for(auto & x : t.vtraj) {
			if (x.first >= endt) break;
			vtraj.insert(x);
		}
	}

	vartraj(vartraj && t) : vtraj(std::move(t.vtraj)) {
		tstart = t.tstart;
		tend = t.tend;
	}

	vartraj & operator =(const vartraj & t) {
		if (this == &t) return *this;
		vtraj = t.vtraj;
		tstart = t.tstart;
		tend = t.tend;
		return *this;
	}

	vartraj & operator =(vartraj && t) {
		vtraj = std::move(t.vtraj);
		tstart = t.tstart;
		tend = t.tend;
		return *this;
	}

	void insert(double t, double v) {
		vtraj.insert(std::make_pair(t, rec(v)));
		if (tstart > t) tstart = t;
		if (tend < t) tend = t;
	}

	void reindex() {}

	double & starttime() {return tstart;}
	const double & starttime() const {return tstart;}
	double & endtime() {return tend;}
	const double & endtime() const {return tend;}

	size_t erase(double key) {return vtraj.erase(key);}

	void print(std::ostream &os) const
	{for(auto & x : *this) os << x.first << ": " << x.second << std::endl;}

private:
	double tstart, tend;
	trajT vtraj;

private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar & ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(tstart);
		ar & BOOST_SERIALIZATION_NVP(tend);
		ar & BOOST_SERIALIZATION_NVP(vtraj);
	}
};
BOOST_CLASS_EXPORT_KEY(vartraj)
BOOST_CLASS_EXPORT_KEY(vartraj::rec)


class traj : public std::vector<vartraj> {
public:
	traj(std::initializer_list<vartraj> data)
	{for(auto & x : data) push_back(x);}
	traj() = default;
	traj(int nvar) : std::vector<vartraj>(nvar) {}
	traj(const traj & tr, double endt) : sx(tr.sx)
	{for(auto & vtr : tr) emplace_back(vtr, endt);}
	std::vector<double> sx;
private:
	typedef std::vector<vartraj> trajectory;
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar & ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(trajectory);
		ar & BOOST_SERIALIZATION_NVP(sx);
	}
};
BOOST_CLASS_EXPORT_KEY(traj)

void printtr(
	std::ostream & os,
	const traj & tr,
	bool incolumns = true,
	bool ishrs = false);

#endif
