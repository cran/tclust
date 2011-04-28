/*
    SMat - Simple Matrix Classes
    Copyright (C) 2011 by Heinrich Fritz (heinrich_fritz@hotmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//	R.meal.h
//	R Mathematical Environment Abstraction Layer

#ifndef R_MEAL_H
#define R_MEAL_H
//#include "R.meal.h"

	class CRmealSettings 
	{
	public:
		CRmealSettings () ;
		CRmealSettings (const char *szEmail) ;
		const char *GetEmail () { return m_szEmail ; }
	protected:
		const char *m_szEmail ;
	} ;

#define R_MEAL_SETTINGS CRmealSettings g_R_meal_settings

#endif	//	#ifndef R_MEAL_H
