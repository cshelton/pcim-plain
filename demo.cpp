#include "pcim.h"
#include <vector>
#include <iostream>
#include <random>
#include <utility>

// a simple program to show how to use the PCIM code

// build the c-cpim from Figure 2 of Parikh, Gunawardana, and Meek
pcim buildmodel() {
	shptr<pcimtest>
		roottst = shptr<pcimtest>(new counttest(1,-1,1,0)),
		ftst = shptr<pcimtest>(new vartest(0)),
		fftst = shptr<pcimtest>(new counttest(1,0,5,0)),
		ttst = shptr<pcimtest>(new counttest(1,-1,2,1));
	
	return pcim(roottst,
			shptr<pcim>(new pcim(ttst,
				shptr<pcim>(new pcim(0.0)),
				shptr<pcim>(new pcim(10.0)))),
			shptr<pcim>(shptr<pcim>(new pcim(ftst,
				shptr<pcim>(new pcim(0.1)),
				shptr<pcim>(new pcim(fftst,
					shptr<pcim>(new pcim(0.1)),
					shptr<pcim>(new pcim(0.0))))
			)))
		);
}

std::vector<traj> sampledata(const pcim &model, int nvar, int ntraj, double T) {
	std::random_device rd;
	std::mt19937 rand(rd());

	std::vector<traj> ret;
	for(int i=0;i<ntraj;i++)
		ret.emplace_back(model.sample(T,nvar,rand));
	return ret;
}

std::vector<shptr<pcimtest>> maketestbank(int nvar) {
	std::vector<shptr<pcimtest>> ret;
	std::vector<int> counts = {1,2,4,8};
	std::vector<std::pair<double,double>> lags = {{1,0},{2,0},{5,0},
										{2,1},{5,1},{5,2}};
	for(int i=-1;i<nvar;i++) {
		for(int c : counts)
			for(auto l : lags)
				ret.emplace_back(new counttest(c,i,l.first,l.second));
		if (i>=0) ret.emplace_back(new vartest(i));
	}
	return ret;
}

pcim learnmodel(const std::vector<traj> &data,
		std::vector<shptr<pcimtest>> testbank) {

	// hyper-parameters set here (chosen pretty arbitrarily
	//  as alpha = 2.5, beta = 2.0, kappa=0.1, minnumevents=5
	//  with the number of threads (numproc) = 4
	pcim::pcimparams hyperparams(25.0,20.0,0.01,50,4);
	return pcim{data,testbank,hyperparams};
}

int main(int argc, char **argv) {

	int nsamp = argc>1 ? atoi(argv[1]) : 10;

	pcim model = buildmodel();
	std::cout << "original model:" << std::endl;
	model.print(std::cout);

	std::vector<traj> data = sampledata(model,2,nsamp,20.0);
	std::cout << "first sampled trajectory:" << std::endl;
	printtr(std::cout,data[0]);

	std::vector<shptr<pcimtest>> testbank = maketestbank(2);
	pcim learnedmodel = learnmodel(data,testbank);


	std::cout << "learned model:" << std::endl;
	learnedmodel.print(std::cout);
}
