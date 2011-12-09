#include "ZedEventsLibrary.h"

ZedEventsLibrary::ZedEventsLibrary() {}

ZedEventsLibrary::~ZedEventsLibrary() {}

ZedEventsLibrary::ZedEventsLibrary(string selection, bool makeFromData = true)
{
    _selection  = selection;
    _dataPeriod = dataPeriod;
    _triggers   = triggers;
}

bool ZedEventsLibrary::GetBin()
{
  return kTRUE;
}
