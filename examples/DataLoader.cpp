#include "DataLoader.h"

namespace analysis {

namespace detail{

  herr_t collectTableNames(hid_t group_id, const char * member_name, void* operator_data){
      std::set<std::string>* items=static_cast<std::set<std::string>*>(operator_data);
      items->insert(member_name);
      return(0);
  }

} // namespace detail

} // namespace analysis
