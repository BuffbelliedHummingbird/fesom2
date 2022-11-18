program open_comfort

	use, intrinsic :: iso_c_binding

	use sqlite3
	use sqlite3_util
	
	implicit none
	
	type(c_ptr)                                :: db
	character(len=40), pointer, dimension(:,:) :: result
	character(len=80)                          :: errmsg
	integer                                    :: rc
	
! Open database
rc = sqlite3_open_v2('file:testDB.db', db, SQLITE_OPEN_READONLY)

! Close database
rc = sqlite3_close(db)


end program

! compile with:
! mpiifort -o open_comfort open_comfort.F90 libfortran-sqlite3.a -lsqlite3
