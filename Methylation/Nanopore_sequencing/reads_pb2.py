# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: reads.proto

import sys
_b=sys.version_info[0]<3 and (lambda x:x) or (lambda x:x.encode('latin1'))
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor.FileDescriptor(
  name='reads.proto',
  package='',
  syntax='proto3',
  serialized_options=None,
  serialized_pb=_b('\n\x0breads.proto\",\n\x04Read\x12\x0c\n\x04uuid\x18\x01 \x01(\t\x12\x16\n\x0emod_base_probs\x18\x02 \x01(\x0c\" \n\x08\x46\x61kReads\x12\x14\n\x05reads\x18\x01 \x03(\x0b\x32\x05.Readb\x06proto3')
)




_READ = _descriptor.Descriptor(
  name='Read',
  full_name='Read',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='uuid', full_name='Read.uuid', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR),
    _descriptor.FieldDescriptor(
      name='mod_base_probs', full_name='Read.mod_base_probs', index=1,
      number=2, type=12, cpp_type=9, label=1,
      has_default_value=False, default_value=_b(""),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=15,
  serialized_end=59,
)


_FAKREADS = _descriptor.Descriptor(
  name='FakReads',
  full_name='FakReads',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='reads', full_name='FakReads.reads', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=61,
  serialized_end=93,
)

_FAKREADS.fields_by_name['reads'].message_type = _READ
DESCRIPTOR.message_types_by_name['Read'] = _READ
DESCRIPTOR.message_types_by_name['FakReads'] = _FAKREADS
_sym_db.RegisterFileDescriptor(DESCRIPTOR)

Read = _reflection.GeneratedProtocolMessageType('Read', (_message.Message,), {
  'DESCRIPTOR' : _READ,
  '__module__' : 'reads_pb2'
  # @@protoc_insertion_point(class_scope:Read)
  })
_sym_db.RegisterMessage(Read)

FakReads = _reflection.GeneratedProtocolMessageType('FakReads', (_message.Message,), {
  'DESCRIPTOR' : _FAKREADS,
  '__module__' : 'reads_pb2'
  # @@protoc_insertion_point(class_scope:FakReads)
  })
_sym_db.RegisterMessage(FakReads)


# @@protoc_insertion_point(module_scope)
