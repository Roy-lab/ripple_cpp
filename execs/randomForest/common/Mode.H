#ifndef _MODE
#define _MODE
class Mode
{
	public:
		Mode();
		~Mode();
		typedef enum
		{
			CLASSIFICATION,
			REGRESSION
		}ModeCode;
		
		static const char* getModeString(int);
	private:
		static const char* modeString[];
};
#endif
