# Class: Inflate

Streaming DEFLATE decompression

## Table of contents

### Constructors

- [constructor](Inflate.md#constructor)

### Properties

- [ondata](Inflate.md#ondata)

### Methods

- [push](Inflate.md#push)

## Constructors

### constructor

• **new Inflate**(`opts`, `cb?`): [`Inflate`](Inflate.md)

Creates a DEFLATE decompression stream

#### Parameters

| Name | Type | Description |
| :------ | :------ | :------ |
| `opts` | [`InflateStreamOptions`](../interfaces/InflateStreamOptions.md) | The decompression options |
| `cb?` | [`FlateStreamHandler`](../README.md#flatestreamhandler) | The callback to call whenever data is inflated |

#### Returns

[`Inflate`](Inflate.md)

• **new Inflate**(`cb?`): [`Inflate`](Inflate.md)

Creates a DEFLATE decompression stream

#### Parameters

| Name | Type | Description |
| :------ | :------ | :------ |
| `cb?` | [`FlateStreamHandler`](../README.md#flatestreamhandler) | The callback to call whenever data is inflated |

#### Returns

[`Inflate`](Inflate.md)

## Properties

### ondata

• **ondata**: [`FlateStreamHandler`](../README.md#flatestreamhandler)

The handler to call whenever data is available

## Methods

### push

▸ **push**(`chunk`, `final?`): `void`

Pushes a chunk to be inflated

#### Parameters

| Name | Type | Description |
| :------ | :------ | :------ |
| `chunk` | `Uint8Array` | The chunk to push |
| `final?` | `boolean` | Whether this is the final chunk |

#### Returns

`void`
