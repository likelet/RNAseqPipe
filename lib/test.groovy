/**
 * Created by likelet on 2019/1/13.
 */


def store = ''
def n=3
n.times{
    store += 'x'
}

static def addstringToalign(String str,int num){
    if(str.length() < num) {
        def numSpace = num - str.length()

        numSpace.times {
            str += ' '
        }
    }
    str

}

def print_parameter(content, parameter){
    println addstringToalign(content, 30)+parameter
}

print_parameter("A", "test")
print_parameter("Aaaaa", "test")
print_parameter("ASSSS   asda", "test")
print_parameter("Addddasdasd  ", "test")