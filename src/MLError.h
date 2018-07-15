#pragma once
#ifndef MUSCLEMASS_SRC_ML_ERROR_H_
#define MUSCLEMASS_SRC_ML_ERROR_H_
#define MLErrorReturn(err) { if (!err.isOK())	{ return err; }}

#include <string>

class MLError
{
public:
	MLError();
	MLError(std::string internalDescription);
	std::string internalDescription() const;
	bool isOK() const;

private:
	bool isOK_;
	std::string internalDescription_;
};

#endif // MUSCLEMASS_SRC_ML_ERROR_H_