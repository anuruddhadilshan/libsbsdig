#include "TSBSSimDetector.h"
#include "TSBSDBManager.h"

TSBSSimDetector::TSBSSimDetector() : fHasData(false)
{
  fDBmanager = TSBSDBManager::GetInstance();
}

TSBSSimDetector::~TSBSSimDetector()
{
}
