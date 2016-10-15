//
//  strtokenizer.h
//  evolda
//
//  Created by huangxin on 15/12/13.
//  Copyright © 2015年 huangxin. All rights reserved.
//

#ifndef strtokenizer_h
#define strtokenizer_h
#include <string>
#include <vector>

using namespace std;

class strtokenizer {
protected:
    vector<string> tokens;
    int idx;
    
public:
    strtokenizer(string str, string seperators = " ");
    
    void parse(string str, string seperators);
    
    int count_tokens();
    string next_token();
    void start_scan();
    
    string token(int i);
};


#endif /* strtokenizer_h */
