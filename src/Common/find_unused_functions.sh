filelines=`egrep "^function " db.php functions.php genomics.php | cut -f2 -d' ' | cut -f1 -d'('`
for line in $filelines ; do
	count=`find /mnt/users/ahsturm1/SVN/ /mnt/users/ahsturm1/Sandbox/ /mnt/share/data/ /mnt/share/chips/ /mnt/share/doc/ /mnt/share/primer/ /mnt/users/all/http/ -name "*.php" | grep -v test_ | grep -v "+old" | xargs grep "$line(" | grep -v "function $line" | wc -l`
	if [ $count -eq 0 ]; then
		echo "HIT $line $count"
	else
		echo "    $line $count"
	fi
done
