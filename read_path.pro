FUNCTION READ_PATH, WRK_PATH = WRK_PATH, SAV_PATH = SAV_PATH, SYS_PATH = SYS_PATH, $
               BP_PATH = BP_PATH, PAR_PATH = PAR_PATH, DATA_PATH = DATA_PATH, MAC = MAC, $
               M3 = M3, TTT = TTT

  If NOT Keyword_set( WRK_PATH ) AND NOT Keyword_set( DATA_PATH ) then DATA_PATH  = 1

  If Keyword_set( WRK_PATH )  then if Keyword_set( M3 ) then PATH = '/wrk/hsergi/' else PATH = '/global/u1/h/hsergi/'
;;; M3/NERSC
;;  If Keyword_set( DATA_PATH ) then PATH = '/data/hsergi/'
;;   If Keyword_set( DATA_PATH ) then PATH = '/space/hsergi/'
  If Keyword_set( DATA_PATH ) then PATH = './'
;  If Keyword_set( DATA_PATH ) then PATH = '/scratch2/scratchdirs/hsergi/'
;;;
  If Keyword_set( TTT ) then PATH = path + 'ttt/'
  If Keyword_set( SAV_PATH ) then  PATH = path + 'save_sets/'
  If Keyword_set( SYS_PATH ) then  PATH = path + 'maps/systematics/'
  If Keyword_set( BP_PATH ) then   PATH = path + 'polar/'
  If Keyword_set( PAR_PATH ) then  PATH = path + 'par/'
  If Keyword_set( MAC ) then PATH = '/Users/hildebrandt/0WRK/HCal/Processing/Labtools/SH/data/'
return, path
END
