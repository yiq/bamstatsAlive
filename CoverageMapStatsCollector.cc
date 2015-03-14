#include "CoverageMapStatsCollector.h"

using namespace BamstatsAlive;
using namespace std;

CoverageMapStatsCollector::CoverageMapStatsCollector(const GenomicRegionStore::GenomicRegionT * currentRegion, const coverageHistT& existingHistogram) : 
	AbstractStatCollector(), 
	_currentRegion(currentRegion), 
	_coveredLength(0), 
	_existingCoverageHist(existingHistogram)
{
	auto regionLength = _currentRegion->endPos - _currentRegion->startPos + 1;
	_regionalCoverageMap = new unsigned int [regionLength];
	memset(_regionalCoverageMap, 0, sizeof(unsigned int) * regionLength);
}

CoverageMapStatsCollector::~CoverageMapStatsCollector() {
	delete [] _regionalCoverageMap;
}

void CoverageMapStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
	if(!_currentRegion->contains(_currentRegion->chrom, al.Position) && !_currentRegion->contains(_currentRegion->chrom, al.Position + al.Length))
		return;

	auto readMappedStartPos = al.Position - _currentRegion->startPos;
	auto readMappedEndPos = al.Position + al.Length - _currentRegion->startPos;
	auto regionLength = _currentRegion->endPos - _currentRegion->startPos + 1;


	if(readMappedStartPos < 0) readMappedStartPos = 0;
	if(readMappedEndPos >= regionLength) readMappedEndPos = regionLength - 1;

	for(auto i=readMappedStartPos; i <= readMappedEndPos; i++) _regionalCoverageMap[i]++;


	if(readMappedStartPos > _coveredLength) {
		for(size_t i=_coveredLength; i<readMappedStartPos; i++) {
			auto cov = _regionalCoverageMap[i];
			auto covEntryIter = _coverageHist.find(cov);
			if(covEntryIter != _coverageHist.end())
				_coverageHist[cov]++;
			else
				_coverageHist[cov] = 1;
		}
		_coveredLength = readMappedStartPos;
	}
}

CoverageMapStatsCollector::coverageHistT CoverageMapStatsCollector::getEffectiveHistogram(unsigned int& totalPos) {
	totalPos = 0;
	coverageHistT effHist;
	for(auto it = _coverageHist.cbegin(); it != _coverageHist.cend(); it++) {
		if(_existingCoverageHist.find(it->first) != _existingCoverageHist.cend())
			effHist[it->first] = _existingCoverageHist.at(it->first) + it->second;
		else
			effHist[it->first] = it->second;

		totalPos += effHist[it->first];
	}

	return effHist;
}

void CoverageMapStatsCollector::appendJsonImpl(json_t * jsonRootObj) {
	// Coverage Histogram
	json_t * j_cov_hist = json_object();

	unsigned int totalPos = 0;
	coverageHistT effHist = getEffectiveHistogram(totalPos);

	for(auto it = effHist.begin(); it != effHist.end(); it++) {
		stringstream labelSS; labelSS << it->first;
		json_object_set_new(j_cov_hist, labelSS.str().c_str(), json_real( it->second / static_cast<double>(totalPos)));
	}
	json_object_set_new(jsonRootObj, "coverage_hist", j_cov_hist);
}

