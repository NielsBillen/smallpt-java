default: all

all: clean jar

clean:
	$(info removing all .class files)
	@rm -f src/*.class
	$(info removing smallpt-java.jar)
	@rm -f smallpt-java.jar

jar:
	$(info compiling src/SmallPT.java)
	@javac -g -d src -classpath src src/SmallPT.java
	$(info compiling smallpt-java.jar)
	@jar -cfe smallpt-java.jar SmallPT -C src/ .
