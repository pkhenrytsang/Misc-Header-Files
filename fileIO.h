/**********************

Scientific data file IO operations (tested with -std=c++17)

Feb 24 2021

by Pak Ki Henry Tsang

Please do not redistribute without permission from author

email: henrytsang222@gmail.com

**********************/

#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include<vector>
#include<tuple>
#include<cassert>
#include<type_traits>
#include<iterator>


namespace fileIO{

	/*
		Utility Functions
	*/

	namespace utility{
		//Read data from file and store into string
		std::string read_data(const char* FN);
		
		//Operation on data extracted from read_data(), one should use these functions to minimize IO (at cost of memory)
		std::vector<std::string> get_comments(std::string data, const char prefix);
		void count_lines(std::string data,size_t* nrows,size_t* ncols,const size_t skiplines);
		
		//Operate on a file directly
		std::vector<std::string> get_comments(const char* FN, const char prefix); //Output a vector of comments (if present)
		void count_lines(const char* FN,size_t* nrows,size_t* ncols,const size_t skiplines); //Counts the number of lines in a file, one should first use get_comments to see if there are comments or not.
		
		template<typename T> std::vector<T> read_col_template(std::string data,const size_t col,const size_t skiplines);
	
		template<class T,class... Targs> std::vector<T> get_row(const size_t idx,Targs const & ... args);
		template<class... Targs> std::vector<size_t> get_size(Targs const & ... args);
		template<typename T,class... Ts> void append(std::vector<T> const& values,Ts & ... args);
		template<class... Ts,class T, size_t... Idxs> void push_back(std::vector<T> const& values, Ts & ... args, std::index_sequence<Idxs...>);
	}
	
	/*
		Column reading
	*/
	//Read all columns (otherwise there will be error)
	template<class... Ts> void read(const char* FN,Ts & ... args);
	template<class... Ts> void read(std::string FN,Ts & ... args){read(FN.c_str(),args...);}
	
	
	//Read specific column
	template<typename T> std::vector<T> read_col(const char* FN,const size_t col);
	template<typename T> std::vector<T> read_col(std::string FN,const size_t col){return read_col<T>(FN.c_str(),col);}
	
	//Read first column
	template<typename T> std::vector<T> read_col(const char* FN);
	template<typename T> std::vector<T> read_col(std::string FN){return read_col<T>(FN.c_str());}
	
	
	/*
		Column writing
	*/
	
	template<class... Targs> void write(const char* FN,Targs const & ... args);
	template<class... Targs> void write(std::string FN,Targs const & ... args){write(FN.c_str(),args...);};
	
};

/*******************************************************************
	Read data into memory - first step before any further operations
*******************************************************************/

std::string fileIO::utility::read_data(const char* FN){

	//Reads data into a string
	//Taken from https://stackoverflow.com/questions/2912520/read-file-contents-into-a-string-in-c
	
	using namespace std;
	
  ifstream ifs(FN);
  string data( (istreambuf_iterator<char>(ifs) ),
                       (istreambuf_iterator<char>()    ) );
                       
  return data;
}

/******************************
	Read Column
*******************************/

//specific col : char* version
template<typename T> std::vector<T> fileIO::read_col(const char* FN,const size_t col){

	using namespace std;
	
	size_t nrows,ncols;
	
	string data = utility::read_data(FN); //Read data
	
	vector<string> comments = utility::get_comments(data, '#'); //Read in comments with default delimiter
	utility::count_lines(data,&nrows,&ncols,comments.size()); //Count number of lines in file
	if (col>ncols){	cout<<"parameter col exceed number of columns in file\n"<<endl ;abort();} //if column to read exceeds that in file, abort
	
	vector<T> output = utility::read_col_template<T>(data,col,comments.size());
	
	return output;
}

//first col : char* version
template<typename T> std::vector<T> fileIO::read_col(const char* FN){

	using namespace std;
	
	size_t nrows,ncols;
	size_t col=1;
	
	string data = utility::read_data(FN); //Read data
	
	vector<string> comments = utility::get_comments(data, '#'); //Read in comments with default delimiter
	utility::count_lines(data,&nrows,&ncols,comments.size()); //Count number of lines in file
	if (col>ncols){	cout<<"parameter col exceed number of columns in file\n"<<endl ;abort();} //if column to read exceeds that in file, abort
	
	vector<T> output = utility::read_col_template<T>(data,col,comments.size());
	
	return output;
}

//Routine template
template<typename T> std::vector<T> fileIO::utility::read_col_template(std::string data,const size_t col,const size_t skiplines){

	using namespace std;
	
	vector<T> output;
	
	istringstream iss_data(data);
	
	string line; //buffer
	size_t i=0;
	while(i<skiplines && iss_data.peek()!=EOF)	{getline(iss_data, line); i++;} //skip lines
	
	size_t row_idx=1;
	//Begin column counting
	while (iss_data.peek() != EOF){
	
		getline(iss_data, line); //read line to buffer
		row_idx++;
		
		size_t col_idx=1;
		
		istringstream iss(line);
		do{
			string dummy;
			iss >> dummy;
			if (dummy.length()) {if (col_idx==col) output.push_back( (T) stod(dummy)); col_idx++;}
		}
		while(iss);
	}
	
	return output;
}

/*
template<class... Ts,class T, size_t... Idxs>
std::tuple<Ts...>
fileIO::parse(std::vector<T> const& values, std::index_sequence<Idxs...>) {
    return {values[Idxs]...};
}

template<class... Ts,class T>
std::tuple<Ts...> fileIO::create_tuple(std::vector<T> const& values) {
    assert(sizeof...(Ts) == values.size());
    return parse<Ts...>(values, std::make_index_sequence<sizeof...(Ts)>{});
}
template<class... Ts>
std::tuple<Ts...> fileIO::readf(const char* FN) {

	using namespace std;
	
	using T = std::common_type_t<Ts...>;
	using F = typename T::value_type;
	
	vector<T> columns;
	size_t nrows,ncols;
	
	string data = read_data(FN); //Read data
	
	vector<string> comments = get_comments(data, '#'); //Read in comments with default delimiter
	count_lines(data,&nrows,&ncols,comments.size()); //Count number of lines in file
	
	istringstream iss_data(data);
	string line; //buffer
	size_t i=0;
	while(i<comments.size() && iss_data.peek()!=EOF)	{getline(iss_data, line); i++;} //skip lines
	
	
	for (size_t col = 1 ;col<=ncols;col++) columns.push_back( T() );
	
	//Save column data into vector of doubles
	//while (iss_data.peek() != EOF){
	for (size_t row = 1; row<=nrows ;row++){
	
		getline(iss_data, line); //read line to buffer

		istringstream iss(line);
		
		for (size_t col = 1; col<=ncols ;col++){ //iterate over columns
			string dummy;
			iss >> dummy;
			if (dummy.length()) columns[col-1].push_back( (F) stod(dummy));
		}
	}

	return create_tuple<Ts...>(columns);
}
*/

template<class... Ts>
void fileIO::read(const char* FN,Ts & ... args) {

	/*
		Reads in variable number of columns
	*/

	using namespace std;
	
	using T = std::common_type_t<Ts...>;
	using F = typename T::value_type;
	
	size_t nrows,ncols;
	
	string data = utility::read_data(FN); //Read data
	
	vector<string> comments = utility::get_comments(data, '#'); //Read in comments with default delimiter
	utility::count_lines(data,&nrows,&ncols,comments.size()); //Count number of lines in file
	
	istringstream iss_data(data);
	string line; //buffer
	size_t i=0;
	while(i<comments.size() && iss_data.peek()!=EOF)	{getline(iss_data, line); i++;} //skip lines
	
	//while (iss_data.peek() != EOF){
	for (size_t row = 1; row<=nrows ;row++){
	
		getline(iss_data, line); //read line to buffer

		istringstream iss(line);
		T row_data;
		for (size_t col = 1; col<=ncols ;col++){ //iterate over columns
			string dummy;
			iss >> dummy;
			if (dummy.length()) row_data.push_back( (F) stod(dummy));
		}
		utility::append(row_data,args...);
	}
}

template<class... Ts,class T, size_t... Idxs>
void fileIO::utility::push_back(std::vector<T> const& values, Ts & ... args, std::index_sequence<Idxs...>) {
	((args.push_back(values[Idxs])), ...);
}

template<typename T,class... Ts>
void fileIO::utility::append(std::vector<T> const& values,Ts & ... args) {
	using namespace std;
  assert(sizeof...(Ts) == values.size());
	push_back<Ts...>(values, args...,make_index_sequence<sizeof...(Ts)>{});
}

/*
	Write to file
*/

template<class T,class... Targs> 
std::vector<T> fileIO::utility::get_row(const size_t idx,Targs const & ... args) {
		std::vector output={args[idx]...};
    return output;
}

template<class... Targs> 
std::vector<size_t> fileIO::utility::get_size(Targs const & ... args) {
		std::vector output={args.size()...};
    return output;
}

template<class... Targs> 
void fileIO::write(const char* FN,Targs const & ... args){

	using namespace std;
	using T = common_type_t<Targs...>;
	using F = typename T::value_type;
	
	auto sizes = utility::get_size(args...);
	assert( equal(sizes.begin() + 1, sizes.end(), sizes.begin()) ); //Check if all elements equal
	
	const size_t nrows=sizes[0];
	
	ofstream file;
	
	file.open(FN);
	if(file.fail()) {
    cout << "error opening file" << endl; abort();}
  {
		for (size_t row=0;row<nrows;row++){
			auto tmp = utility::get_row<F>(row,args...);
			
			ostringstream oss;

			if (!tmp.empty())
			{
				copy(tmp.begin(), tmp.end()-1,ostream_iterator<F>(oss, " "));
				oss << tmp.back();
			}

			file << oss.str() << endl;
		}
	}file.close();

}


/********************
	Utility Functions
*********************/

std::vector<std::string> fileIO::utility::get_comments(const char* FN, const char prefix){

	using namespace std;
	
	vector<string> comments;
	
	ifstream file;
	
	file.open(FN);
	if(file.fail()) {
    cout << "error opening file" << endl; 
    abort();
	}{
		while(file.peek() != EOF){
			string line;
			getline(file, line); //read line to buffer
			if(line.at(0)==prefix){ line.erase(0,1); comments.push_back(line);}
			else{break;}
		}
	}file.close();
	
	
	return comments;
}

void fileIO::utility::count_lines(const char* FN,size_t* nrows,size_t* ncols,const size_t skiplines){

	/*
		If 
	*/

	using namespace std;
	
	ifstream file;
	
	//Set # of rows and cols to zero
	*nrows=0;
	*ncols=0;
	
	//Count # of cols of first row
	
	file.open(FN);
	if(file.fail()) {
    cout << "error opening file" << endl; 
    abort();
	}{
		string line; //buffer
		size_t i=0;
		while(i<skiplines && file.peek()!=EOF)	{getline(file, line); i++;} //skip lines
		
		//Begin column counting
		if (file.peek() != EOF){
			getline(file, line); //read line to buffer
			istringstream iss(line);
			do{
				string dummy;
				iss >> dummy;
				if (dummy.length()) *ncols=*ncols+1;
			}
			while(iss);
		}
		else{ cout << "no data left after skipping lines" << endl;  abort(); }
	}file.close();
	
	//Now we are sure there are data to read, we count the number of rows
	file.open(FN);{
		string dummy; //buffer
		for(size_t i=0;i<skiplines;i++)	getline(file, dummy); //skip lines
	
		while(file.peek() != EOF){
			getline(file, dummy);
			*nrows=*nrows+1;
		}
	
	}file.close();
}

std::vector<std::string> fileIO::utility::get_comments(std::string data, const char prefix){

	using namespace std;
	
	vector<string> comments;
	
	istringstream iss(data);
	
	while(iss.peek() != EOF){
		string line;
		getline(iss, line); //read line to buffer
		if(line.at(0)==prefix){ line.erase(0,1); comments.push_back(line);}
		else{break;}
	}
	
	return comments;
}

void fileIO::utility::count_lines(std::string data,size_t* nrows,size_t* ncols,const size_t skiplines){

	/*
		If 
	*/

	using namespace std;
	
	istringstream iss_col(data);{
		//Set # of rows and cols to zero
		*nrows=0;
		*ncols=0;
		
		//Count # of cols of first row
		string line; //buffer
		size_t i=0;
		while(i<skiplines && iss_col.peek()!=EOF)	{getline(iss_col, line); i++;} //skip lines
		
		//Begin column counting
		if (iss_col.peek() != EOF){
			getline(iss_col, line); //read line to buffer
			istringstream iss(line);
			do{
				string dummy;
				iss >> dummy;
				if (dummy.length()) *ncols=*ncols+1;
			}
			while(iss);
		}
		else{ cout << "no data left after skipping lines" << endl;  abort(); }
	}
	
	istringstream iss_row(data);{
		//Now we are sure there are data to read, we count the number of rows
		string dummy; //buffer
		for(size_t i=0;i<skiplines;i++)	getline(iss_row, dummy); //skip lines

		while(iss_row.peek() != EOF){
			getline(iss_row, dummy);
			*nrows=*nrows+1;
		}
	}
}


