program open_comfort

	use sqlite3
	
	implicit none
	
	type(SQLITE_DATABASE)                      :: db
	character(len=40), pointer, dimension(:,:) :: result
	character(len=80)                          :: errmsg


end program

! compile with:
! mpiifort -o open_comfort open_comfort.F90 libfortran-sqlite3.a -lsqlite3
