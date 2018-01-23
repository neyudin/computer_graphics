all:
	make cross
	
build:
	cd ../template; make all
clean:
	cd ../template; make clean	
train:
	../template/build/bin/task2 -d ./data/multiclass/train_labels.txt -m ./tmp/model.txt --train
predict:
	../template/build/bin/task2 -d ./data/multiclass/test_labels.txt -m ./tmp/model.txt -l ./tmp/predictions.txt -p 
	
	
CUR_PATH = $(shell pwd)

test:
	python ./src/my_compare.py $(CUR_PATH)/data/multiclass/test_labels.txt $(CUR_PATH)/tmp/predictions.txt ./test_output.txt

cross:
	make build
	
	rm -f output.txt
	echo 'tests ' | cat > output.txt
	
	for i in `seq 49`; do \
		python ./src/gendata.py .5 ; \
		make train ; \
		make predict ; \
		python ./src/my_compare.py ./data/multiclass/test_labels.txt ./tmp/predictions.txt ./output.txt; \
	done
	
	python ./src/average.py ./output.txt
