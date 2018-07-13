#include "MLError.h"

MLError::MLError()
{
	isOK_ = true;
	internalDescription_ = "no error";
}
MLError::MLError(std::string internalDescription)
{
	isOK_ = false;
	internalDescription_ = internalDescription;
}
std::string MLError::internalDescription() const
{
	return internalDescription_;
}
bool MLError::isOK() const
{
	return isOK_;
}
