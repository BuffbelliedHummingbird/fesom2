PROGRAM pingpong

	USE mpi
	
	IMPLICIT NONE
	
	INTEGER :: ierror, sizecomm
	INTEGER, allocatable :: mymessage(:)
	INTEGER :: counter, count, dest, tag, rank
	INTEGER :: stats(MPI_STATUS_SIZE)
	REAL :: time1, time2, time3, time4, time5, time6
	INTEGER :: myarray(5)
	
	CALL MPI_INIT(ierror)
	
	open (unit = 10, file = "pingpongtimer.dat")

	
	DO counter = 1, 10
	! count = 1
	count = 10**counter
	allocate(mymessage(count))

	CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERROR)
	
	tag = 0
	
	IF (rank==0) THEN
	
		time1 = MPI_Wtime()
		call MPI_SEND(mymessage, count, MPI_INT, 1, tag, MPI_COMM_WORLD, ierror)
		!time2 = MPI_Wtime()
		
		!WRITE(*,*) 'Sending time (1): ', time2-time1
	END IF
	
	IF (rank==1) THEN
		!time3 = MPI_Wtime()
		call MPI_RECV(mymessage, count, MPI_INT, 0, tag, MPI_COMM_WORLD,stats,ierror)
		!time4 = MPI_Wtime()
		
		call MPI_SEND(mymessage, count, MPI_INT, 0, tag, MPI_COMM_WORLD, ierror)
		
		!WRITE(*,*) 'Receiving time (2): ', time4-time3
	END IF
	
	IF (rank==0) THEN
		call MPI_RECV(mymessage, count, MPI_INT, 1, tag, MPI_COMM_WORLD,stats,ierror)
		time5 = MPI_Wtime()

		WRITE(10,*) time5-time1
	END IF
	
	deallocate(mymessage)
	
	END DO
	close(10)

	myarray(:) = 5
	
	write(*,*) myarray
	
	CALL MPI_FINALIZE(ierror)

END

! compile with:
! mpif90 -o hello hello.F90 

! running the program:
! mpirun -np TASKS hello
! (TASKS is a number specifying the number of processes)
