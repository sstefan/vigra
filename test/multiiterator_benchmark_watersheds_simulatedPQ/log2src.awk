#!/usr/bin/awk -f
BEGIN {
  FS="[[:space:],]+";
}
/^_________BEGIN_________/ {
    print "#include <cstring>"
    print "struct LogData {"
    print "  static const size_t shape[];"
    print "  static const size_t data[];"
    print "  static const size_t count;"
    print "};"
    printf("const size_t LogData::data[] = {\n");
    count=0;
}
/^__SHAPE:/ {
    shape1=$2;
    shape2=$3;
    shape3=$4;
}
/^__push:/ {
    printf(" 1, %d, %d,\n", $2, $3);
    count+=1;
}
/^__pop:/ {
    printf(" 0, %d, %d,\n", $2, $3);
    count+=1;
}
/^_________END_________/ {
    printf("};\n\n");
    printf(" const  size_t LogData::shape[] = {%d, %d, %d};\n",shape1, shape2, shape3);
    printf(" const  size_t LogData::count = %d;\n",count);
}

END{}
