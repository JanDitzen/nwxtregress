*! nwxtregress
*! v. a0.03
*! 30.08.2021 - see https://janditzen.github.io/nwxtregress/

capture program drop nwxtregress
program define nwxtregress, eclass
syntax varlist(ts min=2) [if] 	, 	[	///
		/// spatial weight matrix settings
		dvarlag(string)					/// name of spatial weights matrix for dependent var
		seed(string)					/// set seed	
		* 								/// ivlag weights
		update							/// update from github
		version							/// shows version
		*								/// remaining options
		] 
		

		if "`update'" != "" {
			qui nwxtregress, version
			local v_installed "`r(version)'"	
			cap net uninstall nwxtregress
			if _rc == 0 {
				noi disp "Version `v_installed' removed."
				local rm = 0
			}
			else {
				cap ssc uninstall nwxtregress
				if _rc == 0 {
					noi disp "Version `v_installed' (from SSC) removed."
				}
				else {
					noi disp "No version from net install or ssc install found!"
				}
			}
			noi disp "Updating from Github ... ", _c
			cap qui net install nwxtregress , from("https://janditzen.github.io/nwxtregress/") force replace
			if _rc == 0 {
				noi disp " Update successfull."
				qui nwxtregress, version
				noi disp "New version is `r(version)'"
			}
			else {
				noi disp "Update not successfull!"
			}
			exit
		}

		if "`version'" != "" {
			local version 0.03
			noi disp "This is version `version' - 30.08.2021"
			nwxtregress , version
			return local version "`version'"
			exit
		}

		**** Main Program

		noi disp as text ""
		noi disp as error "This is an alpha version! Results may change.", _c
		noi disp as text ""

		*** necessary checks
		* moremata
		*** check that moremata is installed
		qui{
			cap findfile moremata.hlp
			if _rc != 0 {
				noi disp "Moremata is required. Please install:"
				noi disp in smcl "{stata ssc install moremata}"
				error 198
			}
		}

		* set seed
		if "`seed'" != "" {
			set seed `seed'
		}

		if  regexm("`options'","ivarlag") == 1 & "`dvarlag'" != "" {
			_nwxtregsdm `varlist' `if' , `options' dvarlag(`dvarlag')		
		}
		else if regexm("`options'","ivarlag") == 0 & "`dvarlag'" != "" {
			_nwxtregsar `varlist' `if' ,  `options' dvarlag(`dvarlag')
		}
		else if regexm("`options'","ivarlag") == 0 & "`dvarlag'" == "" {
			/// OLS MODEL - to be included
			noi disp "OLS model not supported yet"
		}
		else {
			noi disp "Option ivarlag() or dvarlag() required."
		}

end
