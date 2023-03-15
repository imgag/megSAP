filelines=`egrep "^function " db.php functions.php genomics.php | cut -f2 -d' ' | cut -f1 -d'('`
for line in $filelines ; do
	count=`find /mnt/storage3/users/ahsturm1/ngs-bits/ /mnt/storage3/users/ahsturm1/megSAP/ /mnt/storage3/users/ahsturm1/Sandbox/ /mnt/storage2/megSAP/data/ /mnt/storage1/share/chips/ /mnt/storage1/share/doc/ /mnt/storage1/share/primer/ /mnt/storage1/users/bioinf/http/ -name "*.php" | grep -v test_ | grep -v "+old" | xargs grep "$line(" | grep -v "function $line" | wc -l`
	if [ $count -eq 0 ]; then
		echo "HIT $line $count"
	else
		echo "    $line $count"
	fi
done
