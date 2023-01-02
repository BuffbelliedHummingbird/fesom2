rm usage.log
touch usage.log
for i in {1..240}
do
	sleep 1
	info.sh -l >> usage.log
done
